% Get a list of the data types in the GaitData structure
% Find the indices of all data types that are not 'meta'
dataType_idx = find(~strcmp(fieldnames(GaitData),'meta'));
% Get a list of GaitData fieldnames
dataTypes = fieldnames(GaitData);
% Filter out the 'meta' fieldname
dataTypes = dataTypes(dataType_idx);


for d = 1:length(dataTypes) % Loop over each data type
    for k = 1:7
        for s = 1:2
            if s == 1
                side = "L";
            else
                side = "R";
            end

            %% BIC
            %k = width(NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W); % n synergies
            Xhat = NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W * NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).H; % Reconstructed Data
            X = GaitData.(dataTypes{d}).(side).concat; % Original activations
            m = length(NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W); % n time points
            n = width(NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).H); % n muscles
            c = min(sqrt(m),sqrt(n));

            BIC1 = log((norm(Xhat-X,'fro')^2)) + k*((m+n)/m*n) * log(m*n/(m+n));
            BIC2 = log((norm(Xhat-X,'fro')^2)) + k*((m+n)/m*n) * log(c^2);
            BIC3 = log((norm(Xhat-X,'fro')^2)) + k*((m+n)/m*n) * (log(c^2))/(c^2);

            %% AIC
            D = GaitData.(dataTypes{d}).(side).concat'; % Original activations transposed (Muscles x Activations)
            C = NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).H'; % Muscle Weightings transposed (Synergies x Muscle Weights)
            W = NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W'; % Activation Patterns transposed (Activations x Synergies)
            N = width(GaitData.(dataTypes{d}).(side).concat); % n Muscles
            M = length(GaitData.(dataTypes{d}).(side).concat); % n time points
            %T = width(NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W); % n synergies
            %AIC = sum((D - (C*W)).^2,'all') + 2*N * (M+T);

            E = sum((D - (C*W)).^2,'all');
            variance = var(GaitData.(dataTypes{d}).(side).concat,0,'all');
            psi = N * (M+k);
            AIC = 2*(E/(2*variance)+psi);

            %AIC = sum((D - (C*W)).^2,'all')/2*var(GaitData.(dataTypes{d}).(side).concat,0,'all') + N * (M+k);

            NNMF.informationCriteria.(dataTypes{d}).(side)(k,1) = BIC1;
            NNMF.informationCriteria.(dataTypes{d}).(side)(k,2) = BIC2;
            NNMF.informationCriteria.(dataTypes{d}).(side)(k,3) = BIC3;
            NNMF.informationCriteria.(dataTypes{d}).(side)(k,4) = AIC;
        end
    end
end