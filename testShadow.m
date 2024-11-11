clear
clc
%%
allRsoFixPositions = getRsoPostions();
sunFixPosArray = getSunFixPositions ( allRsoFixPositions(:,1,4) );
%%
au  = 149597870699.999988; % Astronomical unit [m]; DE440
rs = 696000000.0;
re = 6378137;
angumb = atan((rs-re)/au);
angpen = atan((rs+re)/au);

figure()
rotate3d

[xx,yy,zz] = sphere(500);
mesh(xx*6378137, yy*6378137, zz*6378137)
view(3)
axis([-45168531 45168531 -45168531 45168531 -45168531 45168531])
axis square
grid off
hold on

for i = 1 : size(allRsoFixPositions,1)
    l = line([0 sunFixPosArray(i,1)], [0 sunFixPosArray(i,2)], ...
        [0 sunFixPosArray(i,3)] );
    for j = 1: 76


        reci = [allRsoFixPositions(i,j,1), ...
            allRsoFixPositions(i, j,2), ...
            allRsoFixPositions(i, j,3)];
        rsun = [sunFixPosArray(i,1), ...
            sunFixPosArray(i,2), ...
            sunFixPosArray(i,3)];

        if dot(reci,rsun) < 0.0
            ang1 = angl(-rsun, reci);
            sathoriz = mag(reci )*cos(ang1);
            satvert  = mag(reci )*sin(ang1);
            x = re/sin(angpen);
            penvert = tan(angpen)*(x + sathoriz);
            if satvert <= penvert
                pen = 'y';
                y = re/sin(angumb);
                umbvert = tan(angumb)*(y-sathoriz);
                if satvert <= umbvert
                    umb = 'y';
                    h(j) = plot3(allRsoFixPositions(i,j,1), ...
                        allRsoFixPositions(i, j,2), ...
                        allRsoFixPositions(i, j,3), 'ko');
                else
                    h(j) = plot3(allRsoFixPositions(i,j,1), ...
                        allRsoFixPositions(i, j,2), ...
                        allRsoFixPositions(i, j,3), 'ro');
                end
            else

                h(j) = plot3(allRsoFixPositions(i,j,1), ...
                    allRsoFixPositions(i, j,2), ...
                    allRsoFixPositions(i, j,3), 'ro');


            end


        else
            h(j) = plot3(allRsoFixPositions(i,j,1), ...
                allRsoFixPositions(i, j,2), ...
                allRsoFixPositions(i, j,3), 'ro');

        end


    end




    pause(0.01);
    delete(h);
    delete(l);
end
