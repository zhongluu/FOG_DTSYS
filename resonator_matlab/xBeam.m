classdef xBeam
    % Description:
    % xBeam is consist of various property lights.
    % Nomenclature:
    %   grpLight--a group of lights.
    %   iota--beam intensity. (W / m^2)
    %   symEx--the symbol expression of electric field
    %   symHy--the symbol expression of magnetic field
    % Reference:
    % Ajoy Ghatak.'Optics(6th)'. McGraw-Hill Education and Tsinghua University Press, Chapter 20, pp: 402-404 and
    % Appendix D-E.
    % Declaration:
    % Copyright(c) 2021-2025, by Yulu Zhong, All rights reserved. Southeast University, NanJing, P.R.China
    % 07/31/2021, 07/31/2025
    properties
        grpLight;
        epsilon;
        mu;
        iota;
        symEx;
        symHy;
    end

    properties (Access = private)
        isGausBeam = 0;
        gausOmega0;
        gausLambda;
        gausIota0;
        gausE0;
        gausH0;
    end

    methods
        function obj = xBeam(varargin)
            % lights must be single light or a group of lights
            if nargin == 1
                lights = varargin{1};
                if class(lights) ~= 'xLight'
                    error("xBeam: The input parameter must be a light or an array of light");
                end
                obj.grpLight = lights;
                obj.iota = 0;
                for index = 1 : length(lights)
                    obj.iota = obj.iota + lights(index).getIntensity();
                    obj.symEx = obj.symEx + lights(index).symEx;
                    obj.symHy = obj.symHy + lights(index).symHy;
                end
            else
                switch nargin
                case 5
                    % constructed a Gaussian Beam
                    obj = construcedGaussianBeam(obj, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
                otherwise
                    error("xBeam: Construction input parameters nums error");
                end
            end
        end
        
        function obj = pushLight(obj, lights)
            obj.grpLight = [obj.grpLight, lights];
            for index = 1 : length(lights)
                obj.iota = obj.iota + lights(index).getIntensity();
                obj.symEx = obj.symEx + lights(index).symEx;
                obj.symHy = obj.symHy + lights(index).symHy;
            end
        end

        function obj = deleteLight(obj, indexArr)
            if max(indexArr) > length(obj.grpLight) || min(indexArr) < 1
                error("xBeam: The input index mismatch");
            end
            for index = 1 : length(indexArr)
                obj.iota = obj.iota - obj.grpLight(indexArr(index)).getIntensity();
                obj.symEx = obj.symEx - obj.grpLight(indexArr(index)).symEx;
                obj.symHy = obj.symHy - obj.grpLight(indexArr(index)).symHy;
            end
            obj.grpLight(indexArr) = [];
        end

        function [outEx, outHy] = beaming(obj, inZ, inT)
            if ~(iscolumn(inZ) && iscolumn(inT))
                error("xBeam: Input must be a column vector");
            end
            if length(inZ) ~= length(inT)
                error("xBeam: The length of inZ and inT must be equal");
            end
            outEx = zeros(length(inT), 3, length(obj.grpLight));
            outHy = zeros(length(inT), 3, length(obj.grpLight));
            lights = obj.grpLight;
            parfor index = 1 : length(obj.grpLight)
                [outEx(:, :, index), outHy(:, :, index)] = lights(index).lighting(inZ, inT);
            end
        end

        function outFig = dispMine(obj)
            if ~obj.isGausBeam
                error('xBeam: Please using method dispMineGausBeam() to display Gaussian Beam');
            end
            % pre set -- step: 1000 samples; range: 2 * maxlambda
            minSteps = 1000;
            fullRange = 2;
            % find min and max lambda
            lambda = zeros(1, length(obj.grpLight));
            nu = zeros(1, length(obj.grpLight));
            for index = 1 : length(obj.grpLight)
                lambda(index) = obj.grpLight(index).getWaveLength();
                nu(index) = obj.grpLight(index).getWaveFreq();
            end
            maxLambda = max(lambda);
            minLambda = min(lambda);
            minT = 1.0 / max(nu);
            maxT = 1.0 / min(nu);
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
            Ht = zeros(length(obj.grpLight), length(stepT)); 
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
                Ht(index, :) = tmp.h0 * real(exp(-tmp.omega * stepT * 1i));
                XExz(index, :) = Exz(:, 1)'; YExz(index, :) = Exz(:, 2)'; ZExz(index, :) = Exz(:, 3)';
                XHyz(index, :) = Hyz(:, 1)'; YHyz(index, :) = Hyz(:, 2)'; ZHyz(index, :) = Hyz(:, 3)';
            end
            outFig = figure('Name', 'Beam Property');
            % plot 3D XYZ
            subplot(3, 4, [1, 2, 5, 6]);
            plot3(XExz', YExz', ZExz'); hold on;
            plot3(XHyz', YHyz', ZHyz'); hold off; grid on;
            xlabel('z \it (m)');
            ylabel('Hy \it (A / m)');
            zlabel('Ex \it (V / m)');
            subplot(3, 4, [3, 4, 7, 8]);
            plot3(sum(XExz), sum(YExz), sum(ZExz)); hold on;
            plot3(sum(XHyz), sum(YHyz), sum(ZHyz)); hold off; grid on;
            xlabel('z \it (m)');
            ylabel('Hy \it (A / m)');
            zlabel('Ex \it (V / m)');
            subplot(3, 4, [9, 10]);
            plot(stepT, Et);
            xlabel('t \it (s)');
            ylabel('Et \it (V / m)');
            subplot(3, 4, [11, 12]);
            plot(stepT, Ht); grid on;
            xlabel('t \it (s)');
            ylabel('Ht \it (A / m)');
        end

        function outFig = dispMineGausBeam(obj)
            % H vs t and Hy figure is omit in this method
            if obj.isGausBeam == 0
                error("xBeam: The beam must be a Gaussian beam that construct by xbeam(5params)");
            end
            % pre set -- step: 1000 samples; range: -lambda ~ lambda, -T ~ T
            minSteps = 1000;
            % find min and max lambda
            nu = zeros(1, length(obj.grpLight));
            for index = 1 : length(obj.grpLight)
                nu(index) = obj.grpLight(index).getWaveFreq();
            end
            minT = 1.0 / max(nu);
            maxT = 1.0 / min(nu);
            deltaZ = 6 * obj.gausOmega0 * 10 / minSteps;
            deltaT = minT / minSteps;
            stepT = (-maxT : deltaT : maxT - deltaT)';
            stepZ = (-3 * obj.gausOmega0 : deltaZ : 3 * obj.gausOmega0  - deltaZ)';
            [Exz, ~] = obj.beaming(stepZ, zeros(size(stepZ)));
            outFig = figure('Name', 'Gaussian Beam Property');
            % plot 3D XYZ
            subplot(2, 3, [1, 2, 4, 5]);
            for index = 1 : length(obj.grpLight)
                plot3(Exz(:, 1, index), Exz(:, 2, index), Exz(:, 3, index)); hold on;
            end
            hold off; grid on;
            xlabel('z \it (m)');
            ylabel('Ex \it (V / m)');
            zlabel('Ex \it (V / m)');
            % plot E--t
            [tmpEt, ~] = obj.beaming(zeros(size(stepT)), stepT);
            Et = zeros(length(stepT), length(obj.grpLight));
            for index = 1 : length(obj.grpLight)
                Et(:, index) = sqrt(sum((tmpEt(:, :, index) - obj.grpLight(index).pos).^2, 2));
            end
            subplot(2, 3, 3);
            plot(stepT, Et');
            xlabel('t \it (s)');
            ylabel('Et \it (V / m)');
            % plot the filed distribution
            subplot(2, 3, 6);
            deltaX = 4 * obj.gausOmega0 / minSteps;
            stepX = 2 * obj.gausOmega0  - 4 * obj.gausOmega0 / minSteps : -deltaX : -2 * obj.gausOmega0;
            amplitudeDis = zeros(length(stepX), length(stepZ));
            k = 2 * pi / obj.gausLambda;
            alpha = pi^2 * obj.gausOmega0^4 / obj.gausLambda^2;
            omegaZ = obj.gausOmega0 .* sqrt(1 + stepZ.^2 ./ alpha);
            reRz = stepZ ./ (stepZ.^2 + alpha);
            for i = 1 : length(stepX)
                for j = 1 : length(stepZ)
                    phi = k * stepZ(j) + 0.5 * k * reRz(j) * stepX(i)^2;
                    amplitudeDis(i, j) = real( obj.gausE0 * (obj.gausOmega0 / omegaZ(j)) * exp(-phi * 1i - ...
                                         (stepX(i)^2 / omegaZ(j)^2) + 1i*atan(stepZ(j) / sqrt(alpha))));
                end
            end
            modAmplitudeDis = amplitudeDis.^2;
            imagesc(stepZ, stepX, modAmplitudeDis);
            xlabel('z \it (m)');
            ylabel('x \it (m)');
        end
    end

    methods (Access = private)
        function obj = construcedGaussianBeam(obj, inNumOfR, inNumOfTheta, inOmega0, inLambda, inIota0)
            switch (inNumOfR <= 0) * 2 + (inNumOfTheta <= 0)
            case 0
                numOfR = inNumOfR; numOfTheta = inNumOfTheta;
            case 1
                numOfTheta = 10; numOfR = inNumOfR;
            case 2
                numOfR = 10; numOfTheta = inNumOfTheta;
            otherwise
                % default contain 250 lights
                numOfR = 10; numOfTheta = 10;
            end
            obj.epsilon = 8.8542e-12; % C^2 / (N * m^2)
            obj.mu = 4 * pi * 1e-7; % N * s^2 / C^2
            obj.gausLambda = inLambda;
            obj.gausOmega0 = inOmega0;
            eta =  sqrt(obj.mu / obj.epsilon);
            v = 1 / sqrt(obj.epsilon * obj.mu); % m / s
            k = 2 * pi / inLambda;
            nu = v * k / (2 * pi);
            syms x y z t;
            alpha = pi^2 * inOmega0^4 / inLambda^2;
            omegaz = inOmega0 * sqrt(1 + z^2 / alpha);
            obj.iota = inIota0 / (1 + z^2 / alpha) * exp(-2 * (x^2 + y^2) / omegaz^2);
            % default: 10 samples per circle, 10 samples for the radius, so the whole samples is 10 * 10 = 100
            %          the light range from 3 * omega0 to 0
            obj.grpLight = xLight(numOfR * numOfTheta);
            r = 0 : 3 * inOmega0 / numOfR : 3 * inOmega0 * (numOfR -1) / numOfR;
            theta = 0 : 2 * pi / numOfTheta : 2 * pi * (numOfTheta -1) / numOfTheta;
            for i = 1 : numOfR
                for j = 1 : numOfTheta
                    pos = [0, r(i) * sin(theta(j)), r(i)*cos(theta(j))];
                    polarEx = [0, sin(theta(j)), cos(theta(j))];
                    polarHy = [0, cos(theta(j)), -sin(theta(j))];
                    tmpIota = inIota0 / (1 + z^2 / alpha) * exp(-2 * r(i)^2 / omegaz^2);
                    obj.grpLight((i - 1) * numOfTheta + j) = xLight(tmpIota, nu, pos, polarEx, polarHy);
                end
            end
            % assign syms beam
            Phi = k * z + (k * z / (2 * (z^2 + alpha))) * (x^2 + y^2);
            e0 = sqrt(2 * eta * inIota0);
            h0 = e0 / eta;
            obj.gausE0 = e0;
            obj.gausH0 = h0;
            obj.symEx = [0, 0, 1] * e0 * (obj.gausOmega0 / omegaz) * exp(-(x^2 + y^2) / omegaz) *...
                        exp(1i * (-Phi - v * k * t + atan(z / sqrt(alpha))));
            obj.symHy = [0, 1, 0] * h0 * (obj.gausOmega0 / omegaz) * exp(-(x^2 + y^2) / omegaz) *...
                        exp(1i * (-Phi - v * k * t + atan(z / sqrt(alpha))));
            obj.isGausBeam = 1;
        end
    end

end