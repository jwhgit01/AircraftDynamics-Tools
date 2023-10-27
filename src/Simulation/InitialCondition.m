classdef InitialCondition < handle
    %InitialCondition Summary of this class goes here
    %   Detailed explanation goes here

    properties
        RigidBody % struct
        Unsteady
        Actuator
        Enviroment
        Observer
        Controller
    end

    methods
        function obj = InitialCondition(Simulation)
            %InitialCondition Construct an default instance of this class.

            % Initialize the simulation
            Simulation.Initialize;

            % Rigid Body States
            obj.RigidBody.Position = zeros(3,1);
            obj.RigidBody.Attitude = [];
            obj.RigidBody.Velocity = zeros(3,1);
            obj.RigidBody.AngularVelocity = zeros(3,1);

            % Unsteady States
            obj.Unsteady = zeros(Simulation.Dimensions.n_us,1);

            % Actuator States
            obj.Actuator = zeros(Simulation.Dimensions.n_act,1);

            % Envirmental States
            obj.Enviroment = zeros(Simulation.Dimensions.n_env,1);

            % Observer States
            obj.Observer = zeros(Simulation.Dimensions.n_obsv,1);

            % Controller States
            obj.Controller = zeros(Simulation.Dimensions.n_ctrl,1);

            % Add the IC to the simulation
            Simulation.InitialCondition = obj;

        end

        function x0 = IC(obj)
            x0 = [obj.RigidBody.Position;
                  obj.RigidBody.Attitude;
                  obj.RigidBody.Velocity;
                  obj.RigidBody.AngularVelocity;
                  obj.Unsteady;
                  obj.Actuator;
                  obj.Enviroment;
                  obj.Observer;
                  obj.Controller];
        end
    end

end