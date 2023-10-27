classdef ParameterEstimates < handle & dynamicprops
    %ParameterEstimates Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Name % string
        Mean % vector
        Covariance % matrix
        ParameterNames % cell
        EstimationMethod % string
        TrimCondition % struct
        DatasetName % string (the dataset itself is not very portable)
        ColoredCovariance % Matrix
        Rsquared
        NRMSE
        TIC
    end % properties

    methods
        
        function obj = ParameterEstimates(Name,AeroModel)
            %ParameterEstimates Construct an instance of this class

            % Required object properties
            obj.Name = Name;

            % add to aerodynamic model
            AeroModel.AddParameterEstimates(obj);

        end % ParameterEstimates
        
        function ImportParameterEstimates(obj,FileName)
            %ImportParameterEstimates Import aerodynamic model parameter
            % estimates from CSV file.
            
            % read csv data as table
            opts = detectImportOptions(FileName);
            opts.MissingRule = 'omitvar';
            tbl = readtable(FileName,opts);
            
            % check that data is there
            if height(tbl) < 1
                error('Aerodynamic parameters CSV empty!');
            end

            % check that we at least have parameter estimates
            if ~any(strcmp(tbl.Properties.VariableNames,'Estimate'))
                error('Missing data in parameter estimate column!')
            end
            obj.Mean = tbl.Estimate;

            % add table column if present
            if any(strcmp(tbl.Properties.VariableNames,'Variance'))
                obj.Covariance = diag(tbl.Variance);
            end
            if any(strcmp(tbl.Properties.VariableNames,'ParameterNames'))
                obj.ParameterNames = tbl.ParameterNames;
            end
            if any(strcmp(tbl.Properties.VariableNames,'Regressors'))
                obj.Regressors = tbl.Regressors;
            end     

        end % ImportParameterEstimates

    end % public methods

    methods(Access = protected)
        % Nothing so far
    end % protected methods

    methods(Static)
        
        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = ParameterEstimates(s.Name,s.AeroModel); 
                f = fieldnames(s);
                for ii = 1:length(f)
                    if isprop(newObj,f{ii})
                        newObj.(f{ii}) = s.(f{ii});
                    else
                        addprop(newObj,f{ii});
                        newObj.(f{ii}) = s.(f{ii});
                        mc = metaclass(newObj);
                        warning(['Property ''' f{ii} ''' loaded into object'...
                            'as dynamic property.\nIt should be defined in '...
                            mc.Name '.m or its superclasses\nin order to '...
                            'maintain compatibility.'],[])
                    end
                end
                obj = newObj;
            else
                obj = s;
            end
        end

    end % static methods

end % classdef