function gbObsPos = getGbObsPos()
gbObsPos = zeros(4,3);
Haleakala = [20.07,-156.3,3055];
Nairobi = [-1.17, 36.49, 1795];
Albuquerque = [35.06, -106.36, 1514];
Exmouth = [-21.53, 114.05, 650];
gbObsPos(1,:)=goeToFix(Haleakala);
gbObsPos(2,:)=goeToFix(Nairobi);
gbObsPos(3,:)=goeToFix(Albuquerque);
gbObsPos(4,:)=goeToFix(Exmouth);
end