classdef xLaser
    % Description:
    %xLight is constructed by light intensity and frequency under vaccum condition OR by Maxwell equations with specific
    %electric field intensity and magnetic field intensity.
    % Nomenclature:
    %   e0--electric field intensity.
    %   h0--magnetic field intensity.
    %   omega--angular frequency in plane electromagnetic waves.
    %   ita--medium intrinsic impedance.
    %   k--pattern.
    %   lamtha--light wave length.
    %   ni--light wave frequency.
    %   yiota--light intensity.
    % Reference:
    % Ajoy Ghatak.'Optics(6th)'. McGraw-Hill Education and Tsinghua University Press, Chapter 23, pp: 489-510.
    % Declaration:
    % Copyright(c) 2021-2025, by Yulu Zhong, All rights reserved. Southeast University, NanJing, P.R.China
    % 07/31/2021, 07/31/2025
    
    properties
        Property1
    end
    
    methods
        function obj = xLaser(obj, inNums, inF0Iota, inF0, inRangeNu)
            if inNums <= 0
                % default contain 11 lights
                numOfLights = 11;
            else
                numOfLights = inNums;
            end
            if ~mod(inNums, 2) && inNums ~= 0
                warning("xBeam: Recommend to input a odd num to obtain a best visual");
            end
            sigma = inRangeNu / 2 / 3;
            epsilon = 8.8542e-12; % C^2 / (N * m^2)
            mu = 4 * pi * 1e-7; % N * s^2 / C^2
            v = 1 / sqrt(epsilon * mu);
            stepNu = inF0 - sigma : (sigma / floor(numOfLights / 2)) : inF0 + sigma;
            stepIota = inF0Iota * gaussmf(stepNu, [sigma, inF0]);
            % pre alloc
            obj.grpLight = xLight(length(stepNu));
            for index = 1 : length(stepNu)
                obj.grpLight(index) = xLight(stepIota(index), stepNu(index));
            end
            obj.iota = sum(stepIota);
            syms nu z t;
            gaus = inF0Iota * exp(-(nu - inF0)^2 / (2 * sigma^2));
            obj.symEx = [0, 0, 1] * gaus * exp(1i * (2 * pi * v * z / nu - 2 * pi * nu * t));
            obj.symHy = [0, 1, 0] * gaus * exp(1i * (2 * pi * v * z / nu - 2 * pi * nu * t));
            obj.isGausBeam = 1;
        end
        
        function outputArg = method1(obj)
            % H vs t figure is omit in this method
            if obj.isGausBeam == 0
                error("xBeam: The beam must be a Gaussian beam that construct by xbeam(4params)");
            end
            % pre set -- step: 1000 samples; range: 2 * maxlambda
            minSteps = 1000;
            fullRange = 2;
            % find min and max lambda
            minLambda = obj.grpLight(end).getWaveLength();
            maxLambda = obj.grpLight(1).getWaveLength();
            minT = 1 / obj.grpLight(end).getWaveFreq();
            maxT = 1 / obj.grpLight(1).getWaveFreq();
            deltaZ = minLambda / minSteps;
            deltaT = minT / minSteps;
            stepT = (0 : deltaT : (fullRange * maxT) - deltaT)';
            stepZ = (0 : deltaZ : (fullRange * maxLambda) - deltaZ)';
            % pre alloc
            XExz = zeros(length(obj.grpLight), length(stepZ));
            YExz = zeros(length(obj.grpLight), length(stepZ));
            ZExz = zeros(length(obj.grpLight), length(stepZ));
            XHyz = zeros(length(obj.grpLight), length(stepZ));
            YHyz = zeros(length(obj.grpLight), length(stepZ));
            ZHyz = zeros(length(obj.grpLight), length(stepZ));
            Et = zeros(length(obj.grpLight), length(stepT)); 
            resIota = zeros(1, length(obj.grpLight));
            resNu = zeros(1, length(obj.grpLight));
            for index = 1 : length(obj.grpLight)
                tmp = obj.grpLight(index);
                Z = cross(tmp.X, tmp.Y);            
                Z = repmat(Z, length(stepZ), 1) .* stepZ;
                X = repmat(tmp.X, length(stepZ), 1);
                Y = repmat(tmp.Y, length(stepZ), 1);
                Ex = X .* (tmp.e0 * real(exp(tmp.k * stepZ * 1i)));
                Hy = Y .* (tmp.h0 * real(exp(tmp.k * stepZ * 1i)));
                Exz = Ex + Z;
                Hyz  = Hy + Z;
                Et(index, :) = tmp.e0 * real(exp(-tmp.omega * stepT * 1i));
                XExz(index, :) = Exz(:, 1)'; YExz(index, :) = Exz(:, 2)'; ZExz(index, :) = Exz(:, 3)';
                XHyz(index, :) = Hyz(:, 1)'; YHyz(index, :) = Hyz(:, 2)'; ZHyz(index, :) = Hyz(:, 3)';
                resIota(index) = obj.grpLight(index).getIntensity();
                resNu(index) = obj.grpLight(index).getWaveFreq();
            end
            outFig = figure('Name', 'Gaussian Beam Property');
            % plot 3D XYZ
            subplot(2, 3, [1, 2, 4, 5]);
            plot3(XExz', YExz', ZExz'); hold on;
            plot3(XHyz', YHyz', ZHyz'); hold on;
            % outline the centre frequency
            if ~mod(length(obj.grpLight), 2)
                mid = floor((length(obj.grpLight) + 1) / 2);
                plot3(XExz(mid, :), YExz(mid, :), ZExz(mid, :), 'linewidth', 4); hold on;
                plot3(XHyz(mid, :), YHyz(mid, :), ZHyz(mid, :), 'linewidth', 4); hold on;
                plot3(XExz(mid + 1, :), YExz(mid + 1, :), ZExz(mid + 1, :), 'linewidth', 4); hold on;
                plot3(XHyz(mid + 1, :), YHyz(mid + 1, :), ZHyz(mid + 1, :), 'linewidth', 4); hold on;
            else
                mid = (length(obj.grpLight) + 1) / 2;
                plot3(XExz(mid, :), YExz(mid, :), ZExz(mid, :), 'linewidth', 4); hold on;
                plot3(XHyz(mid, :), YHyz(mid, :), ZHyz(mid, :), 'linewidth', 4); hold on;
            end
            hold off; grid on;
            xlabel('z \it (m)');
            ylabel('Hy \it (A / m)');
            zlabel('Ex \it (V / m)');
            % plot E--t
            subplot(2, 3, 3);
            plot(stepT, Et);
            xlabel('t \it (s)');
            ylabel('Et \it (V / m)');
            % plot iota vs nu Gaussian distribution
            subplot(2, 3, 6);
            plot(resNu, resIota);
            xlabel('\nu \it (Hz)');
            ylabel('\iota \it (W / m^2)');
        end
    end
end

