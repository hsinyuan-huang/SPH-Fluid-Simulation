#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#define EPS 1e-4
#define MX_PARTICLE 1000000

using namespace std;

struct vec{
	float x, y, z;
	vec(){}
	vec(int A){x = (0 == A); y = (1 == A); z = (2 == A);}
	vec(float _x, float _y, float _z){x = _x, y = _y, z = _z;}
	float& operator[](int i){
		if(i == 0) return x;
		if(i == 1) return y;
		if(i == 2) return z;
		return x;
	}
};
vec operator +(vec A, vec B){
	return vec(A.x + B.x, A.y + B.y, A.z + B.z);
}
vec operator -(vec A, vec B){
	return vec(A.x - B.x, A.y - B.y, A.z - B.z);
}
vec operator *(vec A, float c){
	return vec(c * A.x, c * A.y, c * A.z);
}
vec operator *(float c, vec A){
	return vec(c * A.x, c * A.y, c * A.z);
}
vec operator /(vec A, float c){
	return vec(A.x / c, A.y / c, A.z / c);
}
float dot(vec A, vec B){
	return A.x * B.x + A.y * B.y + A.z * B.z;
}
float dist(vec A){
	return sqrt(A.x * A.x + A.y * A.y + A.z * A.z);
}
float dist2(vec A){
	return A.x * A.x + A.y * A.y + A.z * A.z;
}


/* ==== Constants ==== */
const vec g(0, 0, -9.82 / 0.004);

const float mass = 0.00020543;
const float rest_density = 600.0 * 0.004 * 0.004 * 0.004;
const float pdist = pow(mass / rest_density, 1.0/3.0);
const float pradi = 2.0;

const float H = 0.01 / 0.004; // kernel radius
const float H2 = H * H; // h2 = h * h;
const float acc_limit = 20000.0;
const float damping = 256.0;
const float bound_repul = 10000.0;

const float Kp = 5; // Pressure Stiffness
const float visc = 0.2 * 0.004; // Viscosity
const float tension = 2000.0; // Surface Tension
const float dt = 0.0001; // time step

const float Wpoly6C = 315. / 64. / M_PI / pow(H, 9);
const float Grad_WspikeC = 45. / M_PI / pow(H, 6);
const float Lapl_WviscC = 45. / M_PI / pow(H, 6);


/* ==== Useful Functions ==== */
float Wpoly6(vec &r){
	float r2 = dist2(r);
	if(H2 < r2) return 0.f;

	float d = (H2 - r2);
	return Wpoly6C * d * d * d;
}
vec Grad_Wpoly6(vec &r){
	float r2 = dist2(r);
	if(H2 < r2) return vec(0, 0, 0);

	float d = (H2 - r2);
	return Wpoly6C * (-6) * d * d * r;
}
float Lapl_Wpoly6(vec &r){
	float r2 = dist2(r);
	if(H2 < r2) return 0.f;

	return Wpoly6C * (-6) * (H2 - r2) * (3 * H2 - 7 * r2);
}

vec Grad_Wspike(vec &r){
	float r2 = dist2(r);
	if(H2 < r2) return vec(0, 0, 0);

	float rlen = sqrt(r2);
	float d = (H - rlen);
	return (-1) * Grad_WspikeC * d * d / rlen * r;
}

float Lapl_Wvisc(vec &r){
	float rlen = dist(r);
	if(H < rlen) return 0.f;

	return Lapl_WviscC * (H - rlen);
}

float pressure(float density){
	return Kp * (pow(density / rest_density, 7) - 1);
}

int part_n;
vec pos[MX_PARTICLE];
vec vel[MX_PARTICLE];
vec acc[MX_PARTICLE];
float den[MX_PARTICLE];

float bound[3] = {100, 30, 100}; // from 0, 0, 0
vector<int> grid[110][110][110];
float len[3];

void add_particles(vec a1, vec a2){
	float d = pdist;
	for(float x = a1.x + EPS; x <= a2.x - EPS; x += d)
		for(float y = a1.y + EPS; y <= a2.y - EPS; y += d)
			for(float z = a1.z + EPS; z <= a2.z - EPS; z += d){
				pos[part_n] = vec(x, y, z);
				vel[part_n] = vec(0, 0, 0);
				part_n ++;
			}
}

void construct_grid(){
	for(int i = 0; i < part_n; i++){
		int gx = (int)(pos[i][0] / len[0]);
		int gy = (int)(pos[i][1] / len[1]);
		int gz = (int)(pos[i][2] / len[2]);
		grid[gx][gy][gz].push_back(i);
	}
}

void initial_scene(){
	part_n = 0;
	add_particles(vec(0, 0, 0), vec(30, 30, 40));
	construct_grid();

	printf("Particle Number: %d\n", part_n);
}

void calc_density(){
	for(int i = 0; i < part_n; i++){
		int gx = (int)(pos[i][0] / len[0]);
		int gy = (int)(pos[i][1] / len[1]);
		int gz = (int)(pos[i][2] / len[2]);

		den[i] = 0;

		for(int dx = -1; dx <= 1; dx++){
			if(gx + dx < 0) continue;
			for(int dy = -1; dy <= 1; dy++){
				if(gy + dy < 0) continue;
				for(int dz = -1; dz <= 1; dz++){
					if(gz + dz < 0) continue;

					for(auto j : grid[gx + dx][gy + dy][gz + dz]){
						// Assume unit mass
						vec direc = pos[j] - pos[i];
						den[i] += mass * Wpoly6(direc);
					}
				}
			}
		}
	}

}

void obtain_forces(){
	for(int i = 0; i < part_n; i++){
		vec f_out = den[i] * g; // Gravitation

		vec f_tens = vec(0, 0, 0); // Surface Tension
		vec f_pres = vec(0, 0, 0); // Pressure Gradient
		vec f_visc = vec(0, 0, 0); // Viscosity

		int gx = (int)(pos[i][0] / len[0]);
		int gy = (int)(pos[i][1] / len[1]);
		int gz = (int)(pos[i][2] / len[2]);

		for(int dx = -1; dx <= 1; dx++){
			if(gx + dx < 0) continue;
			for(int dy = -1; dy <= 1; dy++){
				if(gy + dy < 0) continue;
				for(int dz = -1; dz <= 1; dz++){
					if(gz + dz < 0) continue;
					for(auto j : grid[gx + dx][gy + dy][gz + dz]){
						// Assume unit mass
						if(j == i) continue; // When j == i, both grad and laplacian vanish

						vec direc = pos[i] - pos[j];

						float ave_P = (pressure(den[i]) + pressure(den[j])) / 2.f;
						f_tens = f_tens - tension * Wpoly6(direc) * den[i] * direc; // Surface Tension
						f_pres = f_pres - ave_P * mass / den[j] * Grad_Wspike(direc); // Pressure Gradient
						f_visc = f_visc + (visc * mass / den[j] * Lapl_Wvisc(direc)) * (vel[j] - vel[i]); // Viscosity
					}
				}
			}
		}

		acc[i] = (f_out + f_pres + f_visc + f_tens) / den[i];
	}

}

void update_position(){
	for(int i = part_n - 1; i >= 0; i--){
		// Empty grid
		int gx = (int)(pos[i][0] / len[0]);
		int gy = (int)(pos[i][1] / len[1]);
		int gz = (int)(pos[i][2] / len[2]);
		grid[gx][gy][gz].pop_back();

		// Clamp back to acceleration limit with direction fixed,
		// to increase stability
		float accel2 = dist2(acc[i]);
		if(accel2 > acc_limit * acc_limit)
			acc[i] = acc[i] / sqrt(accel2) * acc_limit;

		// Bounce back from solid objects (1. reversing velocity)
		/*
		   vel[i] = vel[i] + acc[i] * dt;
		   vec nx = pos[i] + vel[i] * dt;
		   for(int j = 0; j < 3; j++)
		   if(nx[j] < 0 || nx[j] > bound[j])
		   vel[i][j] *= -0.83;
		   pos[i] = pos[i] + vel[i] * dt;
		*/

		// Repulse back from solid objects (2. spring force with damping)
		for(int j = 0; j < 3; j++){
			vec normal = vec(0, 0, 0);
			float xdisp = 0;

			if((pos[i][j] - 0) < pradi){
				normal = vec(j);
				xdisp = (pradi - (pos[i][j] - 0));
			}
			if((bound[j] - pos[i][j]) < pradi){
				normal = vec(j) * -1;
				xdisp = (pradi - (bound[j] - pos[i][j]));
			}
			acc[i] = acc[i] + bound_repul * xdisp * normal
				- damping * dot(vel[i], normal) * normal;
		}

		vel[i] = vel[i] + acc[i] * dt;
		pos[i] = pos[i] + vel[i] * dt;
	}

	construct_grid();
}

FILE *matlab = fopen("matlab_render_wc", "w");
FILE *output = fopen("simulated_wcsph", "w");

int main(int argc, char **argv){
	if(argc != 2){
		fprintf(stderr, "./WCSPHfluid ROUND\n");
		exit(-1);
	}
	int ROUND = atoi(argv[1]);

	for(int i = 0; i < 3; i++)
		len[i] = max(bound[i] / 100.f, H);

	initial_scene();

	fprintf(output, "%d %g %g %g\n", part_n, mass, pradi, H);

	// Integration
	for(int T = 0; T < ROUND; T ++){
		// Output position (Matlab)
		for(int i = 0; i < part_n; i++) if(T % 80 == 0)
			fprintf(matlab, "%f %f %f\n", pos[i].x, pos[i].y, pos[i].z);

		calc_density();

		// Output position and density
		for(int i = 0; i < part_n; i++) if(T % 80 == 0)
			fprintf(output, "%g %g %g %g\n", pos[i].x, pos[i].y, pos[i].z, den[i]);
		
		obtain_forces();
		update_position();
	}
}
