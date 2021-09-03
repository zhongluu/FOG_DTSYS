classdef xLight
    % Description:
    % xLight is constructed by light intensity and frequency under vaccum condition OR by Maxwell equations with 
    % specific electric field intensity and magnetic field intensity. Both symbolic solution and numeric solution are 
    % supported in this object
    % Nomenclature:
    %   e0--electric field intensity. (V / m)
    %   h0--magnetic field intensity. (A / m)
    %   omega--angular frequency in plane electromagnetic waves. (rad / s)
    %   eta--medium intrinsic impedance. (Ohm)
    %   k--pattern. (m)
    %   lambda--light wave length. (m)
    %   nu--light wave frequency. (Hz)
    %   iota--light intensity. (W / m^2)
    %   epsilon--dielectric constant (C^2 / (N * m^2))
    %   mu--magnetic conductivity (N * s^2 / C^2)
    %   x--the unit vector of electric field polarization
    %   y--the unit vector of magnetic field polarization
    %   symEx--the symbol expression of electric field
    %   symHy--the symbol expression of magnetic field
    % Reference:
    % Ajoy Ghatak.'Optics(6th)'. McGraw-Hill Education and Tsinghua University Press, Chapter 23, pp: 489-510.
    % Declaration:
    % Copyright(c) 2021-2025, by Yulu Zhong, All rights reserved. Southeast University, NanJing, P.R.China
    % 07/31/2021, 07/31/2025
    properties
        e0; % V / m
        h0; % A / m
        omega; % rad / s
        k; % m
        epsilon; % C^2 / (N * m^2)
        mu; % N * s^2 / C^2
        X;
        Y;
        symEx;
        symHy;
        pos;
    end

    methods
        function obj = xLight(varargin)
            % default constructed by light intensity and frequency under vaccum condition
            % so the default invoked by xLight(iota, nu);
            obj.epsilon = 8.8542e-12; % C^2 / (N * m^2)
            obj.mu = 4 * pi * 1e-7; % N * s^2 / C^2
            obj.X = [0, 0, 1];
            obj.Y = [0, 1, 0];
            % default position is origin
            obj.pos = [0, 0, 0];
            if nargin == 2
                iota = varargin{1}; nu = varargin{2};
                eta =  sqrt(obj.mu / obj.epsilon);
                obj.e0 = sqrt(2 * eta * iota); % V / m
                obj.h0 = obj.e0 / eta; % A / m
                v = 1 / sqrt(obj.epsilon * obj.mu); % m / s
                lambda = v / nu; % m
                obj.omega = 2 * pi * nu; % rad / s
                obj.k = 2 * pi / lambda; % m
                % build symbol expression
                syms z t;
                obj.symEx =pos + obj.X * obj.e0 * exp(1i * (obj.k * z - obj.omega * t));
                obj.symHy =pos + obj.Y * obj.h0 * exp(1i * (obj.k * z - obj.omega * t));
            else
                switch nargin
                    case 0
                        % for pre alloc
                    case 1
                        % for pre alloc
                        in = varargin{1};
                        if length(in) > 2
                            error("xLight: input specify as [r, c]");
                        end
                        if length(in) == 1
                            m = 1;
                            n = in;
                        else
                            m = in(1);
                            n = in(2);
                        end
                        obj(m, n) = obj;
                    case 3
                        obj = construcedWithPos(obj, varargin{1}, varargin{2}, varargin{3});
                    case 4
                        % constructed by Maxwell
                        obj = construcedByMaxwell(obj, varargin{1}, varargin{2}, varargin{3}, varargin{4});
                    case 5
                        obj = construcedWithPosPolarization(obj, varargin{1}, varargin{2}, varargin{3}, varargin{4},...
                                                            varargin{5});
                    otherwise
                        error("xLight: construction input parameters nums error");
                end
            end
        end

        function outVel = getVelocity(obj)
            outVel = obj.omega / obj.k; % m / s
        end

        function outIota = getIntensity(obj)
            v = 1 / sqrt(obj.epsilon * obj.mu);
            outIota = 0.5 * obj.epsilon * v *  obj.e0^2; % W / m^2
        end   
        
        function outEta = getImpedance(obj)
            outEta =  sqrt(obj.mu / obj.epsilon); % Ohm
        end

        function outLambda = getWaveLength(obj)
            outLambda = 2 * pi / obj.k; % m
        end

        function outNu = getWaveFreq(obj)
            outNu = obj.omega / (2 * pi); % Hz
        end
        
        function obj = setPolarization(obj, inEx, inHy)
            if dot(inEx, inHy) ~= 0
                error("xLight: vector Ex and Hy must be orthogonally");
            end
            obj.X = inEx;
            obj.Y = inHy;
            obj.updateSymbol();
        end

        function obj = setPattern(obj, inK)
            obj.k = inK;
            obj.updateSymbol();
        end

        function [outEx, outHy] = lighting(obj, inZ, inT)
            if ~(iscolumn(inZ) && iscolumn(inT))
                error("xLight: Input must be a column vector");
            end
            if length(inZ) ~= length(inT)
                error("xLight: The length of inZ and inT must be equal");
            end
            Z = cross(obj.X, obj.Y);
            syms z t;
            resZ = repmat(Z, length(inZ), 1) .* inZ;
            outEx = double(real(subs(obj.symEx, {z, t}, {inZ, inT}))) + resZ;
            outHy = double(real(subs(obj.symHy, {z, t}, {inZ, inT}))) + resZ;
        end

        function outFig = dispMine(obj)
            % pre set -- step: 1000 samples; range: -lambda ~ lambda, -T ~ T
            steps = 1000;
            lambda = 2 * pi / obj.k;
            T = (2 * pi) / obj.omega;
            deltaZ = lambda / steps;
            deltaT = T / steps;
            stepZ = (-lambda : deltaZ : (lambda - deltaZ))';
            stepT = (-T : deltaT : (T - deltaT))';
            [Exz, Hyz] = obj.lighting(stepZ, zeros(size(stepZ)));
            [Et, Ht] = obj.lighting(zeros(size(stepT)), stepT);
            Et = sqrt(sum((Et - obj.pos).^2, 2));
            Ht = sqrt(sum((Ht - obj.pos).^2, 2));
            outFig = figure('Name', 'Light Property');
            % plot 3D XYZ
            subplot(2, 3, [1, 2, 4, 5]);
            plot3(Exz(:, 1), Exz(:, 2), Exz(:, 3)); hold on;
            plot3(Hyz(:, 1), Hyz(:, 2), Hyz(:, 3)); hold off; grid on;
            xlabel('z \it (m)');
            ylabel('Hy \it (A / m)');
            zlabel('Ex \it (V / m)');
            % plot E--t
            subplot(2, 3, 3);
            plot(stepT, Et); grid on;
            xlabel('t \it (s)');
            ylabel('Et \it (V / m)');
            % plot H--t
            subplot(2, 3, 6);
            plot(stepT, Ht, 'Color','#D95319'); grid on;
            xlabel('t \it (s)');
            ylabel('Ht \it (A / m)');
        end
    end

    methods (Access = private)
        function obj = construcedWithPos(obj, inIota, inNu, inPos)
            iota = inIota; nu = inNu; obj.pos = inPos;
            eta =  sqrt(obj.mu / obj.epsilon);
            obj.e0 = sqrt(2 * eta * iota); % V / m
            obj.h0 = obj.e0 / eta; % A / m
            v = 1 / sqrt(obj.epsilon * obj.mu); % m / s
            lambda = v / nu; % m
            obj.omega = 2 * pi * nu; % rad / s
            obj.k = 2 * pi / lambda; % m
            % build symbol expression
            syms z t;
            obj.symEx = obj.pos + obj.X * obj.e0 * exp(1i * (obj.k * z - obj.omega * t));
            obj.symHy = obj.pos + obj.Y * obj.h0 * exp(1i * (obj.k * z - obj.omega * t));
        end

        function obj = construcedWithPosPolarization(obj, inIota, inNu, inPos, inEx, inHy)
            iota = inIota; nu = inNu; obj.pos = inPos;
            obj.X = inEx; obj.Y = inHy;
            eta =  sqrt(obj.mu / obj.epsilon);
            obj.e0 = sqrt(2 * eta * iota); % V / m
            obj.h0 = obj.e0 / eta; % A / m
            v = 1 / sqrt(obj.epsilon * obj.mu); % m / s
            lambda = v / nu; % m
            obj.omega = 2 * pi * nu; % rad / s
            obj.k = 2 * pi / lambda; % m
            % build symbol expression
            syms z t;
            obj.symEx = obj.pos + obj.X * obj.e0 * exp(1i * (obj.k * z - obj.omega * t));
            obj.symHy = obj.pos + obj.Y * obj.h0 * exp(1i * (obj.k * z - obj.omega * t));
        end

        function obj = construcedByMaxwell(obj, inE0, inH0, inK, inOmega)
            obj.e0 = inE0;
            obj.h0 = inH0;
            obj.k = inK;
            obj.omega = inOmega;
            syms z t;
            obj.symEx = obj.pos + obj.X * obj.e0 * exp(1i * (obj.k * z - obj.omega * t));
            obj.symHy = obj.pos + obj.Y * obj.h0 * exp(1i * (obj.k * z - obj.omega * t));
        end

        function obj = updateSymbol(obj)
            % update symbol expression
            syms z t;
            obj.symEx = obj.pos + obj.X * obj.e0 * exp(1i * (obj.k * z - obj.omega * t));
            obj.symHy = obj.pos + obj.Y * obj.h0 * exp(1i * (obj.k * z - obj.omega * t));
        end
    end

end