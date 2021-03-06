function B = goodcolormap(map)

switch map
    case 'wr'
        B(1,:)= [ones(1,64) ones(1,64)];
        B(2,:)= 1:-1/127:0;
        B(3,:)= 1:-1/127:0;
    case 'bwr'
        B(1,:)= [0:1/63:1 ones(1,64)];
        B(2,:)= [0:1/63:1 1:-1/63:0];
        B(3,:)= [ones(1,64) 1:-1/63:0];
    case 'mk'
        B(1,:)= 1:-1/127:0;  %[0:1/63:1 ones(1,64)];
        B(2,:)= [zeros(1,64) zeros(1,64)];
        B(3,:)= 1:-1/127:0; % [0:1/63:1 ones(1,64)];
    case 'wmk'
        B(1,:)= [0:1/63:1 ones(1,64)];
        B(2,:)= [zeros(1,64) zeros(1,64)];
        B(3,:)= [0:1/63:1 ones(1,64)];
    case 'bk'
        B(1,:)= [zeros(1,64) zeros(1,64)];
        B(2,:)= [zeros(1,64) zeros(1,64)];
        B(3,:)= 1:-1/127:0;
end