classdef CustomScheduler < hNRScheduler

%%%%%


    properties
        % MovingAvgDataRateWeight Moving average parameter to calculate the average data rate
        MovingAvgDataRateWeight (1, 1) {mustBeNumeric, mustBeNonempty,...
                     mustBeGreaterThanOrEqual(MovingAvgDataRateWeight, 0),...
                     mustBeLessThanOrEqual(MovingAvgDataRateWeight, 1)} = 0.5;

        % UEsServedDataRate Stores DL and UL served data rate for each UE 
        % N-by-2 matrix where 'N' is the number of UEs. For each UE, gNB
        % maintains a context which is used in taking scheduling decisions.
        % There is one row for each UE, indexed by their respective RNTI
        % values. Each row has two columns with following information:
        % Served data rate in DL and served data rate in UL direction.
        % Served data rate is the average data rate achieved by UE till now
        % and serves as an important parameter for doing proportional fair
        % scheduling
        UEsServedDataRate
        ConfiguredGrants
        %configuredSlots = [2 3 7 8 9]
        GrantCounter=1
        receivedConfigurations

    end

    methods
        function obj = CustomScheduler(simParameters)
            %fprintf("from customScheduler\n");
            %hNRSchedulerProportionalFair Construct an instance of this class
            
            % Invoke the super class constructor to initialize the properties
            obj = obj@hNRScheduler(simParameters);
            %obj.configureInitial();
            % Moving average parameter to calculate the average data rate
            if isfield(simParameters, 'MovingAvgDataRateWeight')
                obj.MovingAvgDataRateWeight = simParameters.MovingAvgDataRateWeight;
            end

            obj.UEsServedDataRate = ones(length(obj.UEs), 2);
        end
        function updateGrants(obj)
            disp('Current update:');
            disp(obj.receivedConfigurations);
            UEs=find(obj.receivedConfigurations>0);
            grants=1;
            
%             for i=1:length(UEs)
%               for n=1:obj.receivedConfigurations(UEs(i))
%                   obj.ConfiguredGrants(grants,:)=[UEs(i) obj.receivedConfigurations(i)];
%                   grants=grants+1;
%               end  
%             end
              
              if ~isempty(UEs)
              for n=1:length(UEs)
                  obj.ConfiguredGrants(grants,:)=[UEs(n) obj.receivedConfigurations(UEs(n))];
                  grants=grants+1;
              end
              else 
                  obj.ConfiguredGrants=[];
              end 
            disp('Assigned grants:');
            disp(obj.ConfiguredGrants);
        end

        function configureInitial(obj)
        %fprintf("Configuring grants\n");
        obj.ConfiguredGrants=[
            1 2 10
            2 4 10
            2 6 10
            ];
        %disp(obj.ConfiguredGrants);
         end
        
        function resourceAssignments = runULSchedulerTDD(obj)
            %runULSchedulerTDD Runs the gNB scheduler to assign uplink resources in TDD mode
            % If current slot has DL symbol at the start (Scheduling is
            % only done in DL time), then (i) Scheduler selects the
            % upcoming slots (which contains UL resources) to be scheduled.
            % The criterion used for selecting these slots to be scheduled
            % is: All the slot with UL resources which cannot be scheduled
            % in the next DL slot based on the PUSCH preparation time
            % capability of the UEs. It ensures that the UL resources are
            % scheduled as close as possible to the actual transmission
            % time, respecting the PUSCH preparation time capability of
            % the UEs(ii) The scheduler assigns UL resources of the
            % selected slots among the UEs and returns all the UL
            % assignments done
            %
            % RESOURCEASSIGNMENTS = runULSchedulerTDD(OBJ) runs the
            % scheduler and returns the resource assignments structure.
            %
            % RESOURCEASSIGNMENTS is a cell array of structures that contains the
            % resource assignments information.
            
            resourceAssignments = {};
            % Scheduling is only done in the slot starting with DL symbol
            if find(obj.DLULSlotFormat(obj.CurrDLULSlotIndex+1, 1) == obj.DLType, 1)
                slotsToBeSched = selectULSlotsToBeSched(obj); % Select the set of slots to be scheduled in this UL scheduler run
                numULGrants = 0;
                for i=1:length(slotsToBeSched)
                    allowAssignment=false;
%                     for j = 1:size(obj.ConfiguredGrants,1)
%                        if obj.ConfiguredGrants(j,2)==slotsToBeSched(i)
%                             allowAssignment=false;
%                        end 
%                     end
                    % Schedule each selected slot
                    %step1
                    slotULGrants = scheduleULResourcesSlot(obj, slotsToBeSched(i));
                    resourceAssignments(numULGrants + 1 : numULGrants + length(slotULGrants)) = slotULGrants(:);
                    numULGrants = numULGrants + length(slotULGrants);
                    obj.receivedConfigurations
                if ~isequal(obj.receivedConfigurations,[0 0 0 0 0 0 0 0 0 0])
                    disp('entering')
                     grant = configuredULGrant(obj, slotsToBeSched(i));

                     for n=1: length(grant)
                     resourceAssignments{numULGrants+1}=grant(n);
                     numULGrants=numULGrants+1;
                     resourceAssignments{numULGrants}
                     end
                     
                 end
                end 
                % Update the next to-be-scheduled UL slot. Next UL
                % scheduler run starts assigning resources this slot
                % onwards
                if ~isempty(slotsToBeSched)
                    % If any UL slots are scheduled, set the next
                    % to-be-scheduled UL slot as the next UL slot after
                    % last scheduled UL slot
                    lastSchedULSlot = slotsToBeSched(end);
                    obj.NextULSchedulingSlot = getToBeSchedULSlotNextRun(obj, lastSchedULSlot);
                end
            end
        end


    end
    
    methods(Access = protected)
            function selectedSlots = selectULSlotsToBeSched(obj)
            %selectULSlotsToBeSched Get the set of slots to be scheduled by UL scheduler (for TDD mode)
            % The criterion used here selects all the upcoming slots
            % (including the current one) containing unscheduled UL symbols
            % which must be scheduled now. These slots can be scheduled now
            % but cannot be scheduled in the next slot with DL symbols,
            % based on PUSCH preparation time capability of UEs (It is
            % assumed that all the UEs have same PUSCH preparation
            % capability).
            
            selectedSlots = zeros(obj.NumSlotsFrame, 1);
            % Calculate how far the next DL slot is
            nextDLSlotOffset = 1;
            while nextDLSlotOffset < obj.NumSlotsFrame % Consider only the slots within 10 ms
                slotIndex = mod(obj.CurrDLULSlotIndex + nextDLSlotOffset, obj.NumDLULPatternSlots);
                if find(obj.DLULSlotFormat(slotIndex + 1, :) == obj.DLType, 1)
                    break; % Found a slot with DL symbols
                end
                nextDLSlotOffset = nextDLSlotOffset + 1;
            end
            nextDLSymOffset = (nextDLSlotOffset * 14); % Convert to number of symbols
            
            % Calculate how many slots ahead is the next to-be-scheduled
            % slot
            if obj.CurrSlot <= obj.NextULSchedulingSlot
                % It is in the current frame
                nextULSchedSlotOffset = obj.NextULSchedulingSlot - obj.CurrSlot;
            else
                % It is in the next frame
                nextULSchedSlotOffset = (obj.NumSlotsFrame + obj.NextULSchedulingSlot) - obj.CurrSlot;
            end
            
            % Start evaluating candidate future slots one-by-one, to check
            % if they must be scheduled now, starting from the slot which
            % is 'nextULSchedSlotOffset' slots ahead
            numSlotsSelected = 0;
            while nextULSchedSlotOffset < obj.NumSlotsFrame
                % Get slot index of candidate slot in DL-UL pattern and its
                % format
                slotIdxDLULPattern = mod(obj.CurrDLULSlotIndex + nextULSchedSlotOffset, obj.NumDLULPatternSlots);
                slotFormat = obj.DLULSlotFormat(slotIdxDLULPattern + 1, :);
                
                firstULSym = find(slotFormat == obj.ULType, 1, 'first'); % Check for location of first UL symbol in the candidate slot
                if firstULSym % If slot has any UL symbol
                    nextULSymOffset = (nextULSchedSlotOffset * 14) + firstULSym - 1;
                    if (nextULSymOffset - nextDLSymOffset) < obj.PUSCHPrepSymDur
                        % The UL resources of this candidate slot cannot be
                        % scheduled in the first upcoming slot with DL
                        % symbols. Check if it can be scheduled now. If so,
                        % add it to the list of selected slots
                        if nextULSymOffset >= obj.PUSCHPrepSymDur
                            numSlotsSelected = numSlotsSelected + 1;
                            temp=mod(obj.CurrSlot + nextULSchedSlotOffset, obj.NumSlotsFrame);
                            %allowAssignment=true;
%                             for i = 1:size(obj.ConfiguredGrants,1)
%                                 if obj.ConfiguredGrants(i,2)==temp
%                                     %fprintf("Slot Number %d is configured\n",temp)
%                                     allowAssignment=false;
%                                 end 
%                             end 
                        %if(allowAssignment==true)
                          selectedSlots(numSlotsSelected) = temp;
                        %end
                        end
                    else
                        % Slots which are 'nextULSchedSlotOffset' or more
                        % slots ahead can be scheduled in next slot with DL
                        % symbols as scheduling there will also be able to
                        % give enough PUSCH preparation time for UEs.
                        break;
                    end
                end
                nextULSchedSlotOffset = nextULSchedSlotOffset + 1; % Move to the next slot
            end
            selectedSlots = selectedSlots(1 : numSlotsSelected); % Keep only the selected slots in the array
        end

                
        function eligibleUEs = getEligibleUEs(obj)
           eligibleUEs= [1,2,3,4,5,6,7,8,9,10];
           for n=1: height(obj.ConfiguredGrants)
           eligibleUEs(eligibleUEs == obj.ConfiguredGrants(n,1)) = [];
           end
        end
        
            function uplinkGrants = scheduleULResourcesSlot(obj, slotNum)
            %scheduleULResourcesSlot Schedule UL resources of a slot
            %   UPLINKGRANTS = scheduleULResourcesSlot(OBJ, SLOTNUM)
            %   assigns UL resources of the slot, SLOTNUM. Based on the UL
            %   assignment done, it also updates the UL HARQ process
            %   context.
            %   
            %   SLOTNUM is the slot number in the 10 ms frame whose UL
            %   resources are getting scheduled. For FDD, all the symbols
            %   can be used for UL. For TDD, the UL resources can stretch
            %   the full slot or might just be limited to few symbols in
            %   the slot.
            %
            %   UPLINKGRANTS is a cell array where each cell-element
            %   represents an uplink grant and has following fields:
            %
            %       RNTI                Uplink grant is for this UE
            %
            %       Type                Whether assignment is for new transmission ('newTx'),
            %                           retransmission ('reTx')
            %
            %       HARQId              Selected uplink HARQ process ID
            %
            %       RBGAllocationBitmap Frequency-domain resource assignment. A
            %                           bitmap of resource-block-groups of
            %                           the PUSCH bandwidth. Value 1
            %                           indicates RBG is assigned to the UE
            %
            %       StartSymbol         Start symbol of time-domain resources
            %
            %       NumSymbols          Number of symbols allotted in time-domain
            %
            %       SlotOffset          Slot-offset of PUSCH assignment
            %                           w.r.t the current slot
            %
            %       MCS                 Selected modulation and coding scheme for UE with
            %                           respect to the resource assignment done
            %
            %       NDI                 New data indicator flag
            
            % Calculate offset of the slot to be scheduled, from the current
            % slot
            allowAssignment=true;
            for i = 1:size(obj.ConfiguredGrants,1)
               if obj.ConfiguredGrants(i,2)==slotNum
                    allowAssignment=false;
               end 
            end
            if(allowAssignment==true)
                if slotNum >= obj.CurrSlot
                    slotOffset = slotNum - obj.CurrSlot;
                else
                    slotOffset = (obj.NumSlotsFrame + slotNum) - obj.CurrSlot;
                end

                % Get start UL symbol and number of UL symbols in the slot
                if obj.DuplexMode == 1 % TDD
                    DLULPatternIndex = mod(obj.CurrDLULSlotIndex + slotOffset, obj.NumDLULPatternSlots);
                    slotFormat = obj.DLULSlotFormat(DLULPatternIndex + 1, :);
                    firstULSym = find(slotFormat == obj.ULType, 1, 'first') - 1; % Index of first UL symbol in the slot
                    lastULSym = find(slotFormat == obj.ULType, 1, 'last') - 1; % Index of last UL symbol in the slot
                    numULSym = lastULSym - firstULSym + 1;
                else % FDD
                    % All symbols are UL symbols
                    firstULSym = 0;
                    numULSym = 14;
                end

                if obj.SchedulingType == 0 % Slot based scheduling
                    % Assignments to span all the symbols in the slot
                    uplinkGrants = assignULResourceTTI(obj, slotNum, firstULSym, numULSym);

                else % Symbol based scheduling
                    if numULSym < obj.TTIGranularity
                        uplinkGrants = [];
                        return; % Not enough symbols for minimum TTI granularity
                    end
                    numTTIs = floor(numULSym / obj.TTIGranularity); % UL TTIs in the slot

                    % UL grant array with maximum size to store grants
                    uplinkGrants = cell((ceil(14/obj.TTIGranularity) * length(obj.UEs)), 1);
                    numULGrants = 0;

                    % Schedule all UL TTIs in the slot one-by-one
                    startSym = firstULSym;
                    for i = 1 : numTTIs
                        TTIULGrants = assignULResourceTTI(obj, slotNum, startSym, obj.TTIGranularity);
                        uplinkGrants(numULGrants + 1 : numULGrants + length(TTIULGrants)) = TTIULGrants(:);
                        numULGrants = numULGrants + length(TTIULGrants);
                        startSym = startSym + obj.TTIGranularity;
                    end
                    uplinkGrants = uplinkGrants(1 : numULGrants);
                end
                else
                   %step2
                   startSym=0; 
                   TTIULGrants = assignULResourceTTIConfigured2(obj, slotNum, startSym, obj.TTIGranularity);
                   uplinkGrants( 1 : length(TTIULGrants)) = TTIULGrants(:);
            end
            end

    end

    methods(Access = private)
            function ULGrantsTTI = assignULResourceTTIConfigured2(obj, slotNum, startSym, numSym)
            %assignULResourceTTI Perform the uplink scheduling of a set of contiguous UL symbols representing a TTI, of the specified slot
            % A UE getting retransmission opportunity in the TTI is not
            % eligible for getting resources for new transmission. An
            % uplink assignment can be non-contiguous, scattered over RBGs
            % of the PUSCH bandwidth
            
            grant= configuredULGrant(obj,slotNum);
            combinedRBGAllocationBitmap = combineRBG(obj,grant);
            % Assignment of resources for retransmissions
            [reTxUEs, RBGAllocationBitmap, reTxULGrants] = scheduleRetransmissionsUL(obj, slotNum, startSym, numSym,  combinedRBGAllocationBitmap);
            ULGrantsTTI = reTxULGrants;
            % Assignment of resources for new transmissions, if there
            % are RBGs remaining after retransmissions. UEs which got
            % assigned resources for retransmissions as well as those with
            % no free HARQ process, are not eligible for assignment
            eligibleUEs = getNewTxEligibleUEs(obj, obj.ULType, reTxUEs,'ul');
            if any(~RBGAllocationBitmap) && ~isempty(eligibleUEs) % If any RBG is free in the TTI and there are any eligible UEs
                [~, ~, newTxULGrants] = scheduleNewTxUL(obj, slotNum, eligibleUEs, startSym, numSym, combinedRBGAllocationBitmap);
                ULGrantsTTI = [ULGrantsTTI;newTxULGrants];
            end
            end
        
         function ULGrantsTTI = assignULResourceTTIConfigured(obj, slotNum, startSym, numSym)
                grant= configuredULGrant(obj,slotNum);
                %change 1
                eligibleUEs= getEligibleUEs(obj);
                [~, ~, ULGrantsTTI] = scheduleNewTxUL(obj, slotNum, eligibleUEs, startSym, numSym, grant.RBGAllocationBitmap );
         end 


        function updateUEServedDataRate(obj, linkType, resourceAssignments)
            %updateUEServedDataRate Update UEs' served data rate based on RB assignments
            
            if linkType % Uplink
                mcsTable = obj.MCSTableUL;
                totalRBs = obj.NumPUSCHRBs;
                rbgSize = obj.RBGSizeUL;
                numDMRS = obj.NumPUSCHDMRS;
            else % Downlink
                mcsTable = obj.MCSTableDL;
                totalRBs = obj.NumPDSCHRBs;
                rbgSize = obj.RBGSizeDL;
                numDMRS = obj.NumPDSCHDMRS;
            end
            
            % Store UEs which got grant
            scheduledUEs = zeros(length(obj.UEs), 1);
            % Update served data rate for UEs which got grant
            for i = 1:length(resourceAssignments)
                resourceAssignment = resourceAssignments{i};
                scheduledUEs(i) = resourceAssignment.RNTI;
                averageDataRate = obj.UEsServedDataRate(resourceAssignment.RNTI ,linkType+1);
                mcsInfo = mcsTable(resourceAssignment.MCS + 1, :);
                % Bits-per-symbol is after considering both modulation
                % scheme and coding rate
                bitsPerSym = mcsInfo(3);
                % Number of RBGs assigned to UE
                numRBGs = sum(resourceAssignment.RBGAllocationBitmap(:) == 1);
                if resourceAssignment.RBGAllocationBitmap(end) == 1 && ...
                        (mod(totalRBs, rbgSize) ~=0)
                    % If last RBG is allotted and it does not have same number of RBs as
                    % other RBGs.
                    numRBs = (numRBGs-1)*rbgSize + mod(totalRBs, rbgSize);
                else
                    numRBs = numRBGs * rbgSize;
                end
                achievedTxBits = obj.getResourceBandwidth(bitsPerSym, numRBs, ...
                    resourceAssignment.NumSymbols - numDMRS);
                ttiDuration = (obj.SlotDuration * resourceAssignment.NumSymbols)/14;
                achievedDataRate = (achievedTxBits*1000)/ttiDuration; % bits/sec
                updatedAverageDataRate = ((1-obj.MovingAvgDataRateWeight) * averageDataRate) + ...
                    (obj.MovingAvgDataRateWeight * achievedDataRate);
                obj.UEsServedDataRate(resourceAssignment.RNTI, linkType+1) = updatedAverageDataRate;
                obj.UEsServedDataRate
            end
            scheduledUEs = nonzeros(scheduledUEs);
            unScheduledUEs = setdiff(obj.UEs, scheduledUEs);
            
            % Update (decrease) served data rate for each unscheduled UE
            for i=1:length(unScheduledUEs)
                averageDataRate = obj.UEsServedDataRate(unScheduledUEs(i) ,linkType+1);
                updatedAverageDataRate = (1-obj.MovingAvgDataRateWeight) * averageDataRate;
                obj.UEsServedDataRate(unScheduledUEs(i), linkType+1) = updatedAverageDataRate;
            end
        end

                
        function  configureGrants(obj,RNTI,slotNum)
            temp=[RNTI slotNum 10]
            obj.ConfiguredGrants=[obj.ConfiguredGrants;temp];
        end
        
        function  releaseGrant(obj,RNTI,slotNum)
            for i=1:length(obj.ConfiguredGrants)
               if obj.ConfiguredGrants(i,1)==RNTI && obj.ConfiguredGrants(i,2)==slotNum
                   y(i)=false;
               else
                   y(i)=true;
               end
                
            end
            obj.ConfiguredGrants=obj.ConfiguredGrants(y,:)
        end
        
        function combinedRBG =combineRBG(obj,grant)
            
            for n=1:length(grant(1).RBGAllocationBitmap)
              if length(grant)>1
                  if grant(1).RBGAllocationBitmap(n) ==1 || grant(2).RBGAllocationBitmap(n)==1
                      combinedRBG(n)=1;
                  else 
                      combinedRBG(n)=0;
                  end 
              else 
                  combinedRBG=grant(1).RBGAllocationBitmap;
              end 
            end
        end
%         function grant =configuredULGrant(obj,slotNum)
%             UE= obj.ConfiguredGrants(find(obj.ConfiguredGrants(:,2)==slotNum),1);
%             %UE=3
%             grant = struct();
%             if (slotNum==7)
%                 grant.RNTI = UE;
%                 grant.Type = 'newTx';
%                 grant.HARQId = 0;
%                 grant.RBGAllocationBitmap = [1 1 0 0 1 1 0 1 0 1 1 0 0];
%                 grant.StartSymbol = 0;
%                 grant.NumSymbols = 4;
%                 grant.SlotOffset = 1;
%                 grant.MCS = 24;
%                 grant.RV=0;
%                 grant.NDI = 1;
%             elseif (slotNum==8 )
%                 grant.RNTI = UE;
%                 grant.Type = 'newTx';
%                 grant.HARQId = 0;
%                 grant.RBGAllocationBitmap = [1 1 0 1 1 1 0 1 1 0 0 0 0];
%                 grant.StartSymbol = 0;
%                 grant.NumSymbols = 4;
%                 grant.SlotOffset = 2;
%                 grant.MCS = 14;
%                 grant.RV=0;
%                 grant.NDI = 1;
%                elseif (slotNum==3)
%                 grant.RNTI = UE;
%                 grant.Type = 'newTx';
%                 grant.HARQId = 0;
%                 grant.RBGAllocationBitmap = [1 1 0 0 1 1 0 1 0 1 1 0 0];
%                 grant.StartSymbol = 0;
%                 grant.NumSymbols = 14;
%                 grant.SlotOffset = 2;
%                 grant.MCS = 13;
%                 grant.RV=0;
%                 grant.NDI = 1;
%                 elseif (slotNum==2 )
%                 grant.RNTI = UE;
%                 grant.Type = 'newTx';
%                 grant.HARQId = 0;
%                 grant.RBGAllocationBitmap = [1 1 0 1 1 1 0 1 1 0 0 0 0];
%                 grant.StartSymbol = 0;
%                 grant.NumSymbols = 4;
%                 grant.SlotOffset = 1;
%                 grant.MCS = 14;
%                 grant.RV=0;
%                 grant.NDI = 1;
%              elseif (slotNum==9 )
%                 grant.RNTI = UE;
%                 grant.Type = 'newTx';
%                 grant.HARQId = 0;
%                 grant.RBGAllocationBitmap = [1 1 0 0 1 1 0 1 0 1 1 0 0];
%                 grant.StartSymbol = 0;
%                 grant.NumSymbols = 4;
%                 grant.SlotOffset = obj.NextULSchedulingSlot - obj.CurrSlot;
%                 grant.MCS = 24;
%                 grant.RV=0;
%                 grant.NDI = 1;
% 
%             end
%         end
        function grant =configuredULGrant(obj,slot)
            %grant = struct();
            for n=1: height(obj.ConfiguredGrants)
            obj.ConfiguredGrants
            PRBs=  obj.ConfiguredGrants(n,2);
            UE= obj.ConfiguredGrants(find(obj.ConfiguredGrants(:,2)==PRBs),1);
            %UE=3
            
            grant(n).RNTI = UE;
            grant(n).Type = 'newTx';
            grant(n).HARQId = 0;
            grant(n).StartSymbol = 0;
            grant(n).NumSymbols = 8;
            grant(n).SlotOffset = 1;
            grant(n).MCS = 24;
            grant(n).RV=0;
            grant(n).NDI = 1;
            if (PRBs==1)
                grant(n).RBGAllocationBitmap = [1 0 0 0 0 0 0 0 0 0 0 0 0];
            elseif (PRBs==2)
                grant(n).RBGAllocationBitmap = [0 1 1 0 0 0 0 0 0 0 0 0 0];
            elseif (PRBs==3)
                grant(n).RBGAllocationBitmap = [0 0 0 1 1 1 0 0 0 0 0 0 0];
                %grant(n).RBGAllocationBitmap = [0 0 0 0 0 0 0 0 0 0 1 1 1];
            elseif (PRBs==4)
                grant(n).RBGAllocationBitmap = [0 0 0 0 0 0 1 1 1 1 0 0 0];
                %grant(n).RBGAllocationBitmap = [0 0 0 0 0 0 0 0 0 0 0 0 0];
            elseif (PRBs==5)
                grant(n).RBGAllocationBitmap = [0 0 0 0 0 0 1 1 1 1 0 0 0]; 
            elseif (PRBs==6)
                grant(n).RBGAllocationBitmap = [0 0 0 1 1 1 1 1 1 0 0 0 0];
            elseif (PRBs==7)
                grant(n).RBGAllocationBitmap = [0 1 1 1 1 1 1 1 0 0 0 0 0];
            elseif (PRBs==8)
                grant(n).RBGAllocationBitmap = [1 1 1 1 1 1 1 1 0 0 0 0 0];
            end
            end
        end 
end
end 