classdef xResonatorCavity
    % Description:
    % the resonator cavity is a polarization-maintaining fiber(PMF) resonator cavity
    % Nomenclature:
    % Reference:
    % Ajoy Ghatak.'Optics(6th)'. McGraw-Hill Education and Tsinghua University Press, Chapter 23, pp: 489-510.
    % Declaration:
    % Copyright(c) 2021-2025, by Yulu Zhong, All rights reserved. Southeast University, NanJing, P.R.China
    % 07/31/2021, 07/31/2025
    
    properties
        Property1
    end
    
    methods
        function obj = xResonatorCavity(inputArg1,inputArg2)
            %RESONATOR_CARVITY 构造此类的实例
            %   此处显示详细说明
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end
end

