classdef hRLModule < handle
    properties
        %Defining the needed parameters
        
        %Number of UEs in the simulation
        numberOfUEs = 10 
        
        %Maximum number of grants per UE
        maxGPerUser=8;
        
        %Maximum numbr of total grants
        maxGrants=8;
        
        %The probability of choosing to explore, exploits most of the time with a small chance of exploring.
        epsilon = 1
        
        %Length of the simulation from the initial-state until the terminal-state
        total_episodes = 10000
        
        %Number of epochs in each episode
        max_steps = 50
        
        %The learning rate
        alpha = 0.1
        
        %exploartion decreasing decay for exponential decreasing
        exploration_decreasing_decay = 0.001

        %minimum of exploration proba
        min_exploration_proba = 0.01
        
        lambda=0.1
        
        %Discount rate
        gamma = 0.90
        
        %All possible actions 
        action_space
        
        %All possible states 
        state_space
        
        %The size of the state matrix
        state_size 
        
        %The size of the state matrix
        action_size
        
        %Value matrix
        Q
        
        %Eligibility traces matrix
        E
        
        %current currentAction/ grants
        currentAction
        
        %Holds the the next action to be taken
        nextAction
        
        %Holds the current state of the system
        currentState
        
        %Holds the next state of the system
        nextState
        
        %Intermediate variable that holds the state
        state
        
        %An array of the current buffer statuses of all 4 users 
        bufferStatuses
        
        %Holds a list of all taken actions
        allActions
        
        %Holds the number of the current steps in an episodes
        stepOfepisode
        
        currentepisode
        
        prevBufferStatuse
        
        wastedResources
        
        listOfActions
        
        rewardsOverTime
        
        Elapsed_time
        
        episodeRewards
        
        final
    end 
    
    methods
        function obj = hRLModule ()
            clc;
            if exist('MeanrewardsForBufferStatus.txt', 'file')==2
            delete('MeanrewardsForBufferStatus.txt');
            end
            if exist('rewardsForBufferStatus.txt', 'file')==2
            delete('rewardsForBufferStatus.txt');
            end
            obj.prevBufferStatuse= [0 0 0 0 0 0 0 0 0 0];
            obj.bufferStatuses= [0 0 0 0 0 0 0 0 0 0];
            %Initialize the state space 
            [obj.state_space, obj.state_size] = initializeStates(obj); 
            %Initialize the action space 
            [obj.action_space,obj.action_size] = initializeActions(obj);
            %Initialize the Q matrix with zeros
            obj.Q = initializeValues(obj); 
            %Initialize the E matrix with zeros
            obj.E = zeros(obj.state_size, obj.action_size); 
            %Assume no uesers get grants in the first run of the simulation
            obj.currentAction = [0 0 0 0 0 0 0 0 0 0]; 

            % run the simulation
            for episode=1: obj.total_episodes
                runSimulation(obj);
                if ~(isempty(obj.state))
                for stepOfepisode=1: obj.max_steps
                    obj.currentState= obj.state;
                    currentStateIndex = find(sum(obj.state_space==obj.currentState,2)==obj.numberOfUEs);
                    obj.currentAction=obj.chooseAction(currentStateIndex);
                    currentActionIndex = find(sum(obj.action_space==obj.currentAction,2)==obj.numberOfUEs);
                    x=['Current State', num2str(obj.currentState,[1 1 1 1])];
                    disp(x)
                     %if (~isequal(obj.currentState,[1 1 1 1])) && (obj.stepOfepisode<=obj.max_steps) 
                       if stepOfepisode==1 && episode==1
                           tic
                       end 
                        x=['step: ',num2str(stepOfepisode), ' of episode:',num2str(episode)];
                       disp(x);
                       obj.nextState = obj.state;
                       nextStateIndex = find(sum(obj.state_space==obj.nextState,2)==obj.numberOfUEs); 
                       reward= obj.getReward();
                       obj.nextAction= obj.chooseAction(nextStateIndex);
                       nextActionIndex = find(sum(obj.action_space==obj.nextAction,2)==obj.numberOfUEs);
                       %Calling the update function, and passing S,A,R,S',A'
                       obj.update(currentStateIndex,nextStateIndex,currentActionIndex,nextActionIndex,reward);
                       currentStateIndex = find(sum(obj.state_space==obj.currentState,2)==obj.numberOfUEs);
                       currentActionIndex = find(sum(obj.action_space==obj.currentAction,2)==obj.numberOfUEs);
                       obj.stepOfepisode=obj.stepOfepisode+1;

                    if stepOfepisode==obj.max_steps && episode==1
                        obj.Elapsed_time = toc
                    end 
                    if stepOfepisode==obj.max_steps
                        obj.final(end+1)= mean(obj.episodeRewards(1,:));
                        obj.episodeRewards=[];
                        obj.currentepisode=obj.currentepisode+1;

                    end

                    if obj.currentepisode==obj.total_episodes

                        obj.final
                        obj.Elapsed_time
                        obj.listOfActions
                        obj


                    end
                     %We update the exploration proba using exponential decay formula 
                     obj.epsilon = obj.min_exploration_proba + (obj.epsilon - obj.min_exploration_proba) * exp(-1.0 * episode * obj.exploration_decreasing_decay)
                     runSimulation(obj);

                end 
                end
            end
        end
        
        function runRL(obj)
            %SARSA algorithm:
            close all;
            obj.currentState= obj.state;
            currentStateIndex = find(sum(obj.state_space==obj.currentState,2)==obj.numberOfUEs);
            obj.currentAction=obj.chooseAction(currentStateIndex);
            currentActionIndex = find(sum(obj.action_space==obj.currentAction,2)==obj.numberOfUEs);
            x=['Current State', num2str(obj.currentState,[1 1 1 1])];
            disp(x)
%             if (~isequal(obj.currentState,[1 1 1 1])) && (obj.stepOfepisode<=obj.max_steps) 
             if obj.stepOfepisode<=obj.max_steps
               if obj.stepOfepisode==1 && obj.currentepisode==2
                   tic
               end 
               x=['step: ',num2str(obj.stepOfepisode), ' of episode:',num2str(obj.currentepisode)];
               disp(x);
               obj.nextState = obj.state;
               nextStateIndex = find(sum(obj.state_space==obj.nextState,2)==obj.numberOfUEs); 
               reward= obj.getReward();
               obj.nextAction= obj.chooseAction(nextStateIndex);
               nextActionIndex = find(sum(obj.action_space==obj.nextAction,2)==obj.numberOfUEs);
               %Calling the update function, and passing S,A,R,S',A'
               obj.update(currentStateIndex,nextStateIndex,currentActionIndex,nextActionIndex,reward);
               currentStateIndex = find(sum(obj.state_space==obj.currentState,2)==obj.numberOfUEs);
               currentActionIndex = find(sum(obj.action_space==obj.currentAction,2)==obj.numberOfUEs);
               obj.stepOfepisode=obj.stepOfepisode+1;
               
             end 
            if obj.stepOfepisode==obj.max_steps && obj.currentepisode==2
                obj.Elapsed_time = toc
            end 
            if obj.stepOfepisode==obj.max_steps && obj.currentepisode<= obj.total_episodes
                obj.final(end+1)= mean(obj.episodeRewards(1,:));
                fid = fopen('MeanrewardsForBufferStatus.txt', 'a+');
                fprintf(fid, '%d\n', mean(obj.episodeRewards(1,:)));
                fclose(fid);
                obj.episodeRewards=[];
                obj.stepOfepisode=1;
                obj.currentepisode=obj.currentepisode+1;
               
            end
            
            if obj.currentepisode==obj.total_episodes
                
                obj.final
                obj.Elapsed_time
                obj.listOfActions
                obj
                
               
            end
            if obj.stepOfepisode<=obj.max_steps && obj.currentepisode<= obj.total_episodes
             %We update the exploration proba using exponential decay formula 
             %change this obj.epsilon = max(obj.min_exploration_proba, exp(-obj.exploration_decreasing_decay*obj.currentepisode));
             runSimulation(obj);
            end
        end 
        
        
%         function update(obj,currentStateIndex,nextStateIndex,currentActionIndex,nextActionIndex,reward) 
%            %Function responsible for calculating and updating the state-action value matrix
%            %as well as the elgibility traces. 
%            QNextstate = obj.Q(nextStateIndex,nextActionIndex);
%            obj.E(currentStateIndex,currentActionIndex) = obj.E(currentStateIndex,currentActionIndex)+1; 
%            delta = reward(1)+(obj.gamma*QNextstate) - obj.Q(currentStateIndex,currentActionIndex); %error delta
%            obj.Q = obj.Q + obj.alpha*delta.*obj.E; %update Q for all pairs
%            obj.E = obj.gamma*obj.lambda*obj.E; %update E for all pairs
%            
%            obj.currentState=obj.nextState;
%            obj.currentAction=obj.nextAction;
%            obj.listOfActions(end +1,: )= obj.currentAction
%            
%            
%         end

 function update(obj,currentStateIndex,nextStateIndex,currentActionIndex,nextActionIndex,reward) 
     
     obj.Q(currentStateIndex, currentActionIndex) = (1-obj.alpha) * obj.Q(currentStateIndex, currentActionIndex) +obj.alpha *(reward + obj.gamma*max(obj.Q(nextStateIndex,:)));
 end 
        
        function Q_Values =initializeValues(obj)
            Q_Values = zeros(obj.state_size, obj.action_size)+1; 
           
        end 
        
        function chosenAction = chooseAction(obj, state)
            %The main function for selecting the next configuration/
            %Action. The method uses epsilon greedy for the agent to both
            %explore and exploit.
            randN = rand(1); %generating random double between 0 and 1
            if(randN > obj.epsilon) %greedy
                row=obj.Q(state,:);
                [maxQ,index] = max(row(row<1)); %get next pos
                chosenAction=obj.action_space(index,:)
                if isempty(chosenAction)
                rng('shuffle');
                index = randi([1 obj.action_size]);
                chosenAction = obj.action_space(index,:);
                end 
            else
                rng('shuffle');
                index = randi([1 obj.action_size]);
                chosenAction = obj.action_space(index,:);
            end
            obj.allActions( end+1, 1:10) = chosenAction;
            obj.allActions( end,11) = find(sum(obj.action_space==chosenAction,2)==obj.numberOfUEs);
            
        end 
        
        
        function reward = getReward(obj)
            %Reward function. The reward is calculated using the values of the buffer
            %statuses as well as the std of these values. 
            new=sum(obj.prevBufferStatuse)-sum(obj.bufferStatuses);
            obj.bufferStatuses
            obj.prevBufferStatuse
            
%             reward = new-(5*obj.wastedResources);%- sum(obj.nextState);
            %reward = -sum(obj.bufferStatuses)-(5*obj.wastedResources)-(std(obj.bufferStatuses)*5);%- sum(obj.nextState);
            reward = -sum(obj.bufferStatuses);
            obj.rewardsOverTime(1,end+1)=reward;
            obj.episodeRewards(end+1)=reward;

            obj.rewardsOverTime(2,end)=5*obj.wastedResources;
            obj.rewardsOverTime(3,end)= -sum(obj.bufferStatuses);
            wasted=obj.wastedResources;
            fid = fopen('rewardsForBufferStatus.txt', 'a+');
            fprintf(fid, '%d\n', reward);
            fclose(fid);
            reward
            wasted
            new
        end 
        
        function obj = runSimulation(obj)
        % A function that envokes the live script/ the simulation to run.
           NRTDDSymbolBasedSchedulingPerformanceEvaluationExample
        end 
  
        function receiveSimulationResults(obj,currentState, wastedResources)
         %This function is called by the simulation at the end of each
         %step. It receives the buffer statuses of all users after running
         %with the new configuration
         obj.prevBufferStatuse
         currentState
          obj.prevBufferStatuse=obj.bufferStatuses;
          obj.bufferStatuses=currentState;
          obj.state = obj.convertState(currentState);
          obj.wastedResources=wastedResources;
          %obj.runRL()
          close all;

        end 
        
        
        function  [combs,sizeS] = initializeStates(obj)
        %Function responsible for initializing the state space, which is
        %all possible states for the buffer statuses of all users. The size
        %of the state space = 3^(Number of users)
        vectors = { [1 2 3], [1 2 3], [1 2 3] , [1,2,3],[1 2 3],[1 2 3],[1 2 3],[1 2 3],[1 2 3],[1 2 3]}; 
        n = numel(vectors); 
        combs = cell(1,n); 
        [combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); 
        combs = cat(n+1, combs{:}); 
        combs = reshape(combs,[],n);
        sizeS = length(combs);
        
        end
        
        function  [Action,sizeA] = initializeActions(obj)
        %Function responsible for initializing the state space, which is
        %all possible states for the buffer statuses of all users. The size
        %of the state space = 3^(Number of users)
        a=[0 1 2 3 4 5];
        vectors = { a, a, a , a,a,a,a,a,a,a}; 
        n = numel(vectors); 
        Action = cell(1,n); 
        [Action{end:-1:1}] = ndgrid(vectors{end:-1:1}); 
        Action = cat(n+1, Action{:}); 
        Action = reshape(Action,[],n);
          NumberOfUncofig=8;
            Action=Action(sum(Action,2)<=obj.maxGrants,:);
            Action=Action(sum((Action==0),2)>=NumberOfUncofig,:);
            b=sort(Action,2);
            Action(any(diff(b(:,[9,10]),[],2)==0,2),:)=[];
%                for i=1:obj.numberOfUEs-1
%                    for t=i+1:obj.numberOfUEs-2
%                Action(Action(:,i)==Action(:,t), :) = [];
%                end
%                end
               
            sizeA = length(Action);
        end 
        
        
%         function [Action, sizeA] = initializeActions(obj,UEs,maxGPerUser)
% 
%            %Function responsible for initializing the action space. The action space 
%            %consists of all possible actions given the limit in the maximum number
%            %of grants per user (obj.maxGPerUser) and the maximum total
%            %number of grants (obj.maxGrants)
%             for a=1:(maxGPerUser+1)
%                 aStart=(a-1)*(maxGPerUser+1)^3+1;
%                 aDuration=(maxGPerUser+1)^3-1;
%                 aEnd=aStart+aDuration;
%                 Action(aStart:aEnd,1)=a-1;
%                 for b=1:(maxGPerUser+1)
%                     bStart=aStart+(b-1)*(maxGPerUser+1)^2;
%                     bDuration=(maxGPerUser+1)^2-1;
%                     bEnd=bStart+bDuration;
%                     Action(bStart:bEnd,2)=b-1;
%                     for c=1:(maxGPerUser+1)
%                         cStart=bStart+(c-1)*(maxGPerUser+1)^1;
%                         cDuration=(maxGPerUser+1)^1-1;
%                         cEnd=cStart+cDuration;
%                         Action(cStart:cEnd,3)=c-1;
%                         for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,4)=d-1;
%                         end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,5)=d-1;
%                          end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,6)=d-1;
%                          end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,7)=d-1;
%                          end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,8)=d-1;
%                          end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,9)=d-1;
%                          end
%                          for d=1:(maxGPerUser+1)
%                             dStart=cStart+(d-1)*(maxGPerUser+1)^0;
%                             dDuration=(maxGPerUser+1)^0-1;
%                             dEnd=dStart+dDuration;
%                             Action(dStart:dEnd,10)=d-1;
%                         end
%                     end
%                 end
%             end
%             NumberOfUncofig=2;
%             Action=Action(sum(Action,2)<=obj.maxGrants,:);
%             Action=Action(sum((Action==0),2)>=NumberOfUncofig,:);
%             b=sort(Action,2);
%             Action(any(diff(b(:,[3,4]),[],2)==0,2),:)=[];
%             sizeA = length(Action);
%             
%         end 


        
        function currentState= convertState(obj,currentState)
            %Function used to discretize the buffer statuses of the users to
            % (1, 2 ,3) => (Green, Yellow, Red)/ The discretization is based on the predicted 
            %behaviors of the users 
             for i = 1 :length(currentState)
              if currentState(i)<0.3
                  currentState(i)=1;
              elseif currentState(i)>= 0.3 &&  currentState(i)<0.9
                  currentState(i)=2;
              elseif currentState(i)>=0.9
                  currentState(i)=3;
              end 
                  
            end
           currentState = currentState.';
        end 
    end 
    
    methods (Static)
       function currentAction = getLatestConfigurations(obj)
           %Static function used by the simulation to get the latest
           %configuration/ grants. This function is called at the start of
           %each frame.
           currentAction= obj.currentAction;
        end 
    end 
end 
