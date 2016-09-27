fid = fopen('matlab_render');
%fid = fopen('matlab_render_wc');
A = fscanf(fid,'%f');
[elem, ~] = size(A);

writerObj = VideoWriter('water_elastic.avi');
open(writerObj);

N = 9261;
T = round(elem / (3 * N));
A = reshape(A, 3, N * T);

for i = 1:N:N*T
    plot3(A(1, i:i+N-1), A(2, i:i+N-1), A(3, i:i+N-1), 'co');
    axis([0 120 0 100 0 100]);
    view(180, 0);
    grid on
    drawnow
    
    frame = getframe;
    writeVideo(writerObj, frame);
end

close(writerObj);