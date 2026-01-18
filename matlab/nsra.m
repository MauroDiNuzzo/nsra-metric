% nsra.m
function score = nsra(meas_deltas,pred_deltas,Options)

    % Scalable Null-Stratified Rank Accuracy (NSRA)

    arguments (Input)
        meas_deltas (:,1) double % Measured signed deltas vs control
        pred_deltas (:,1) double % Predicted signed deltas vs control
        Options.Epsilon (1,1) double = 0 % Soft sign tolerance
        Options.TieScore (1,1) double = 0.5 % Score for ties or ambiguous comparisons
        Options.Method (1,1) string ...
            {mustBeMember(Options.Method,["full-pairwise","rank-reduced"])} = "rank-reduced"
    end

    arguments (Output)
        score (1,1) double % NSRA score
    end

    epsilon = Options.Epsilon;
    tie_score = Options.TieScore;
    method = Options.Method;

    G = numel(meas_deltas);

    % Measured strata
    s = zeros(G,1);
    s(meas_deltas >  epsilon) =  1; % U
    s(meas_deltas < -epsilon) = -1; % D

    if (method == "full-pairwise")

        % Full-pairwise (Reference definition)

        % - Explicit evaluation of all admissible gene pairs
        % - Pairwise correctness defined directly
        % - Computationally quadratic: O(G^2)
    
        % Pairwise measured sign matrix (upper triangle only) 
        % sign_gt[i,j] = sign(s[i] - s[j])
        % but 0 if both are NULL
        diff_gt = s - s.';
        sign_gt = sign(diff_gt);
    
        % Mask NULL-NULL pairs
        null_null = (s == 0) & (s.' == 0);
        sign_gt(null_null) = 0;
    
        % Keep upper triangle only
        mask_ut = triu(true(G),1);
        sign_gt = sign_gt(mask_ut);
    
        % Predicted ordering
        diff_pred = pred_deltas-pred_deltas.';
        sign_pred = sign(diff_pred);
        sign_pred = sign_pred(mask_ut);
    
        % Valid comparisons
        valid = sign_gt ~= 0;
        if ~any(valid)
            score = NaN;
            return
        end
    
        sg = sign_gt(valid);
        sp = sign_pred(valid);
    
        % Scoring
        correct = (sp == sg);
        tie = (sp == 0);
    
        score = (sum(correct)+tie_score*sum(tie))/numel(sg);        

    else

        % Rank-reduced

        % - Exploits order structure induced by signed deltas
        % - Reduces pairwise consistency to rank comparisons
        % - Computationally sub-quadratic: O(G log G) time, O(G) memory
    
        % Counts
        nU = sum(s ==  1);
        nN = sum(s ==  0);
        nD = sum(s == -1);
    
        % Total number of valid comparisons
        total_pairs = nU*(nN+nD)+nN*nD;
        if total_pairs == 0
            score = NaN;
            return
        end
    
        % Sort by predicted deltas
        [p_sorted,idx] = sort(pred_deltas,"ascend");
        s_sorted = s(idx);
    
        % Cumulative counters
        seen_D = 0;
        seen_N = 0;

        correct = 0;
    
        i = 1;
        while i <= G

            % Handle ties in predicted deltas
            j = i;
            while j <= G && p_sorted(j) == p_sorted(i)
                j = j+1;
            end
    
            block = s_sorted(i:j-1);
    
            % U vs (D, N) seen so far
            nU_b = sum(block ==  1);
            nN_b = sum(block ==  0);

            % NULL vs D seen so far
            nD_b = sum(block == -1);
    
            % U vs D
            correct = correct+nU_b*seen_D;
            % U vs N
            correct = correct+nU_b*seen_N;    
            % NULL vs D
            correct = correct+nN_b*seen_D;
    
            % Ties within block
            correct = correct+tie_score*(nU_b*nD_b+nU_b*nN_b+nN_b*nD_b);
    
            % Update counters
            seen_D = seen_D+nD_b;
            seen_N = seen_N+nN_b;
    
            i = j;
        end
    
        score = correct/total_pairs;

    end

end
