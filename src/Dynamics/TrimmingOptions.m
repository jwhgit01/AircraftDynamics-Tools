classdef TrimmingOptions
    %TrimmingOptions Summary of this class goes here
    %   Detailed explanation goes here

    properties
        AerodynamicModel
        ParameterEstimates
        SteadyMotion
        InitialGuess
    end

    properties (Access = protected)
        EquilibriumStates
    end

    methods
        function obj = TrimmingOptions(SteadyMotion)
            %TrimmingOptions Construct a default instance of this class

            % Depending on the desired SteadyMotion, set default initial
            % guess.
            obj.SteadyMotion = SteadyMotion;
            if strcmp(SteadyMotion,'WingsLevelConstantAltitude')
                obj.InitialGuess.Heading_deg = 0;
                obj.InitialGuess.PitchAngle_deg = 0;
                obj.InitialGuess.AngleOfAttack_deg = 0;
                obj.InitialGuess.Airspeed_m_s = 20;
            elseif strcmp(SteadyMotion,'WingsLevelClimb')
                obj.InitialGuess.Heading_deg = 0;
                obj.InitialGuess.PitchAngle_deg = 0;
                obj.InitialGuess.AngleOfAttack_deg = 0;
                obj.InitialGuess.Airspeed_m_s = 20;
                obj.InitialGuess.FlightPathAngle_deg = 0;
            elseif strcmp(SteadyMotion,'ConstantAltitudeTurn')
                obj.InitialGuess.BankAngle_deg = 30;
                obj.InitialGuess.AngleOfAttack_deg = 0;
                obj.InitialGuess.Airspeed_m_s = 20;
                obj.InitialGuess.TurnRate_deg_s = 10;
            elseif strcmp(SteadyMotion,'ClimbingTurn')
            elseif strcmp(SteadyMotion,'SteadyPullUp')
            else
                error('Invalid Steady Motion')
            end
        end

        function list = GetEquilibriumStates(obj)

        end
    end
end