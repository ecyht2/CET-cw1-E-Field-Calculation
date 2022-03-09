%Clearing Everything
clear all;  %Clear Variables
clc;        %Clear Screen
 
% Define Variable
f = 2.4e9;      %Frequency (fixed)
c = 3e8;        %Speed of Light in a Vacumm
beta = 2*pi*f/c;  %Wave number

% Algorithms
% Solving for Line of Sight
tx = 0; ty = 4; tz = 2;                                 % Transmitter Coords in terms of [x, y, z]
rx = 10; ry = 3; rz = 2;                                % Reciever Coords in terms of [x, y, z]
rxv = rx:0.01:22;                                       % Reciever moving from p to q
d1 = distance(tx, ty, tz, rxv, ry, rz);                 % Distance
Ed = e_field(d1, beta);                                 % Direct Rays (line-of-sight)
Edf = dB(Ed);                                           % Direct Rays in terms of dB

% Test
tx1 = 0; ty1 = 8; tz1 = 2;
m = (ry - ty1)./(rxv-tx1);
b = ry - (m.*rxv);
py = 6;
px = (py - b)./m;
pz = 2;
d1 = distance(tx, ty, tz, px, py, pz);
d2 = distance(px, py, pz, rxv, ry, rz);
d3 = distance(tx, ty, tz, rxv, ry, rz);
d4 = distance(tx1, ty1, tz1, rxv, ry, rz);

% Solving for top reflection
N = 5;
topRef = zeros(N, width(rxv));
for n = 1:N
   topRef(n,:) =  e_field_x_order_reflection(tx, ty, tz, rxv, ry, rz, beta, [0, 6; 0, 0], 1, 3, n);
end

% Solving for bottom reflection
botRef = zeros(N, width(rxv));
for n = 1:N
   botRef(n,:) =  e_field_x_order_reflection(tx, ty, tz, rxv, ry, rz, beta, [0, 0; 0, 6], 1, 3, n);
end

% Total E-field
Etotal = dB(Ed + sum(topRef, 1) + sum(botRef, 1));

% Plotting Results
% Line of Sight
plot(rxv, Edf, rxv, dB(topRef));
% Setting Up Labels
xlabel('Distance from Point p to Point q (m)')
ylabel('Electric Field (db)')
% Setting Up Legend
legend('Line of Sight', 'Top 1st Order', 'Top 2nd Order', 'Top 3rd Order', 'Top 4th Order', 'Top 5th Order')

function d = distance(x1, y1, z1, x2, y2, z2)
    % This function find the distance between 2 points
    d = sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2);
end

function dB = dB(x)
    % This function convers the e-field into dB
    dB = 20*log10(abs(x));
end

function EField = e_field(d, beta)
    % This function finds the e-field of a ray without considering the
    % reflection coeficient
    EField = (1./d).*exp(-1j*beta*d);
end

function refCofficient = reflection_coeficient(thethai, thethat, epsilon1, epsilon2)
    % This function finds the reflection coeficient given thethai, thethat
    % and both epsilon
    refCofficient = (cos(thethai) - sqrt(epsilon2/epsilon1)*cos(thethat))./(cos(thethai) + sqrt(epsilon2/epsilon1)*cos(thethat));
end

function thethai = thetha_i(d1, d2, d3)
    % This function finds the angle of incidence of the reflection given by
    % the 3 length of the triangle
    thethai = acos( (d2.^2 + d3.^2 - d1.^2)./(2*d2.*d3) )./2;
end

function thethat = thetha_t(thethai, epsilon1, epsilon2)
    % This function finds the angle of refraction using snell's law
    thethat = asin(sin(thethai)*sqrt(epsilon1)/sqrt(epsilon2));
end

function [ ix, iy ] = intersection_point(m1, b1, m2, b2)
    % This function finds the intersection point between 2 lines
    ix = (b1 - b2)./(m2 - m1);  % x value of intersection
    iy = m1.*ix + b1;           % y value of intersection
end

function [ xRef, yRef, zRef ] = reflect_point(x, y, z, m, b)
    % This function takes in a point x, y, z and reflect it over the line
    % mx + b
    if m == 0                                                           % Horizontal Line
        diff = b - y;                                                   % Finding diffrence in y
        xRef = x;                                                       % x reflected point
        yRef = b + diff;                                                % y reflected point
        zRef = z;                                                       % z reflected point
    elseif m == inf                                                     % Vertical Line
        diff = b - x;                                                   % Finding difference in x
        xRef = b + diff;                                                % x reflected point
        yRef = y;                                                       % y reflected point
        zRef = z;                                                       % z reflected point
    else                                                                % None Horizontal or Vertical Line
        mNorm = -1/m;                                                   % Gradient of normal line
        bNorm = y - mNorm*x;                                            % y-intercept of normal line
        [ ix, iy ] = intersection_point(m, b, mNorm, bNorm); iz = z;    % Finding the intersection point of the line and the normal line
        xDiff = ix - x; yDiff = iy - y; zDiff = iz - z;                 % Finding the x, y and z difference between the point and the intersection
        xRef = ix + xDiff; yRef = iy + yDiff; zRef = iz +  zDiff;       % Fiding the reflected point
    end
end

function [ m, b ] = points_to_line(x1, y1, x2, y2)
    % This function finds the line that intersects both the point
    m = (y2 - y1)./(x2 - x1);   % Gradient of the line
    b = y1 - m.*x1;             % y-intercept of the line
end

function EField = e_field_x_order_reflection(tx, ty, tz, rx, ry, rz, beta, walls, epsilon1, epsilon2, nOrder)
    % This function calculates the total e-field of the ray that went
    % through nOrder of reflection

    % Caculate all the reflected tx
    reflectedPoints = zeros(nOrder, 3);                                                                                             % Creating Initial Matrix
    [ reflectedPoints(1, 1), reflectedPoints(1, 2), reflectedPoints(1, 3) ] = reflect_point(tx, ty, tz, walls(1,1), walls(1,2));    % Populating the first reflection
    for n = 1:nOrder-1                                                                                                              % Looping over nOrders
        [ reflectedPoints(n+1, 1), reflectedPoints(n+1, 2), reflectedPoints(n+1, 3) ] = reflect_point( ...
            reflectedPoints(n, 1), reflectedPoints(n, 2), reflectedPoints(n, 3), ...
            walls(mod(n+1, height(walls)) + 1, 1), walls(mod(n, height(walls)) + 1, 2) ...
            );                                                                                                                      % Adding reflected points
    end

    coordSize = [ numel(tx), numel(ty), numel(tz), numel(rx), numel(ry), numel(rz) ];   % Getting the width of each coords array
    nPoints = max(coordSize);                                                           % Finding the coords array that is the biggest
    % Intersection Points and rx and tx points
    intersectionPointsX = zeros(nOrder+2, nPoints);                                                                               % Creating Initial X coords matrix
    intersectionPointsY = zeros(nOrder+2, nPoints);                                                                               % Creating Initial Y coords matrix
    intersectionPointsZ = zeros(nOrder+2, nPoints);                                                                               % Creating Initial Z coords matrix
    intersectionPointsX(1,:) = rx; intersectionPointsY(1,:) = ry; intersectionPointsZ(1,:) = rz;                                    % Reciever Point
    intersectionPointsX((nOrder + 2), :) = tx; intersectionPointsY((nOrder + 2),:) = ty; intersectionPointsZ((nOrder + 2),:) = tz;  % Transmitter Point
    for n = 2:nOrder+1                                                                                                                                          % Looping over all possible intersection points
        [ m, b ] = points_to_line( ...
            intersectionPointsX(n-1,:), intersectionPointsY(n-1), ...
            reflectedPoints(nOrder - n + 2, 1), reflectedPoints(nOrder - n + 2, 2) ...
            );  % Finding Line that connects the image and the targeted point
        cWallID = mod(nOrder - n - 1, height(walls)) + 1;                                                                                                       % Choosing which wall the line intersects with
        [ intersectionPointsX(n,:), intersectionPointsY(n,:) ] = intersection_point(m, b, walls(cWallID, 1), walls(cWallID, 2)); intersectionPointsZ(n, :) = 2; % Saving the intersection point
    end

    % Calculating incident angle and angle of refraction
    thethai = zeros(nOrder, nPoints);         % Creating thethai matrix
    thethat = zeros(nOrder, nPoints);         % Creating thethat matrix
    refCoefficient = zeros(nOrder, nPoints);  % Creating reflection coefficient matrix
    for n = 1:nOrder                                                                                                                                                                        % Looping over nOrder
        d1 = distance(intersectionPointsX(n,:), intersectionPointsY(n,:), intersectionPointsZ(n,:), intersectionPointsX(n+2,:), intersectionPointsY(n+2,:), intersectionPointsZ(n+2,:));    % Calculating edge for the angle
        d2 = distance(intersectionPointsX(n+1,:), intersectionPointsY(n+1,:), intersectionPointsZ(n+1,:), intersectionPointsX(n+2,:), intersectionPointsY(n+2,:), intersectionPointsZ(3,:));% Finding the distance of one of the edge
        d3 = distance(intersectionPointsX(n,:), intersectionPointsY(n,:), intersectionPointsZ(n,:), intersectionPointsX(n+1,:), intersectionPointsY(n+1,:), intersectionPointsZ(n+1,:));    % Finding the distance of one of the edge
        thethai(n,:) = thetha_i(d1, d2, d3);                                                                                                                                                % Finding thethai
        thethat(n,:) = thetha_t(thethai(n,:), epsilon1, epsilon2);                                                                                                                          % Finding thethat
        refCoefficient(n,:) = reflection_coeficient(thethai(n,:), thethat(n,:), epsilon1, epsilon2);                                                                                        % Finding the reflection coefficient                                                             % Updating the shared edge
    end

    d = distance(rx, ry, rz, reflectedPoints(nOrder, 1), reflectedPoints(nOrder, 2), reflectedPoints(nOrder, 3));   % Distance the ray traveled in total
    totalRefCoefficient = prod(refCoefficient, 1);                                                                  % Product of all refflection coefficient
    EField = e_field(d, beta).*totalRefCoefficient;                                                                 % Calculating e-field
end