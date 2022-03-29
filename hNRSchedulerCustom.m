classdef hNRSchedulerCustom < hNRScheduler
    %hNRSchedulerProportionalFair Implements proportional fair scheduler

    %   Copyright 2020 The MathWorks, Inc.

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
        
        % Change 
        receivedConfigurations
        
        % Change 
        ConfiguredGrants
        
        % Change 
        prev =1;
        
        counter =1;
        configurations;
    end

    methods
        
       % Change 
       function updateGrants(obj)
            UEs=find(obj.receivedConfigurations>0);
            grants=1;
              if ~isempty(UEs)
              for n=1:length(UEs)
                  obj.ConfiguredGrants(grants,:)=[UEs(n) obj.receivedConfigurations(UEs(n))];
                  grants=grants+1;
                  obj.configurations=obj.ConfiguredGrants;
              end
              else 
                  obj.ConfiguredGrants=[];
              end 
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

                    %step1
                    slotULGrants = scheduleULResourcesSlot(obj, slotsToBeSched(i));
                    resourceAssignments(numULGrants + 1 : numULGrants + length(slotULGrants)) = slotULGrants(:);
                    numULGrants = numULGrants + length(slotULGrants);
                if ~isequal(obj.receivedConfigurations,[0 0 0 0 0 0 0 0 0 0])
                     grant = configuredULGrant(obj, slotsToBeSched(i));
                    
                     for n=1: length(grant)
                     resourceAssignments{numULGrants+1}=grant(n);
                     numULGrants=numULGrants+1;
                     %resourceAssignments{numULGrants};
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
            % Change
%             for i=1 : length(resourceAssignments)
%                 resourceAssignments{i}
%             end  
            
        end

        
        function obj = hNRSchedulerCustom(simParameters)
            %hNRSchedulerProportionalFair Construct an instance of this class

            % Invoke the super class constructor to initialize the properties
            obj = obj@hNRScheduler(simParameters);

            % Moving average parameter to calculate the average data rate
            if isfield(simParameters, 'MovingAvgDataRateWeight')
                obj.MovingAvgDataRateWeight = simParameters.MovingAvgDataRateWeight;
            end

            obj.UEsServedDataRate = ones(length(obj.UEs), 2);
        end
        

        function [selectedUE, mcsIndex] = runSchedulingStrategy(obj, schedulerInput)
            %runSchedulingStrategy Implements the proportional fair scheduling
            %
            %   [SELECTEDUE, MCSINDEX] = runSchedulingStrategy(OBJ, SCHEDULERINPUT) runs
            %   the proportional fair algorithm and returns the UE (among the eligible
            %   ones) which wins this particular resource block group, along with the
            %   suitable MCS index based on the channel conditions. This function gets
            %   called for selecting a UE for each RBG to be used for new transmission
            %   i.e. once for each of the remaining RBGs after assignment for
            %   retransmissions is completed. According to PF scheduling
            %   strategy, the UE which has maximum value for the PF weightage, i.e. the
            %   ratio: (RBG-Achievable-Data-Rate/Historical-data-rate), gets the RBG.
            %
            %   SCHEDULERINPUT structure contains the following fields which scheduler
            %   would use (not necessarily all the information) for selecting the UE to
            %   which RBG would be assigned.
            %
            %       eligibleUEs    - RNTI of the eligible UEs contending for the RBG
            %       RBGIndex       - RBG index in the slot which is getting scheduled
            %       slotNum        - Slot number in the frame whose RBG is getting scheduled
            %       RBGSize        - RBG Size in terms of number of RBs
            %       cqiRBG         - Uplink Channel quality on RBG for UEs. This is a
            %                        N-by-P  matrix with uplink CQI values for UEs on
            %                        different RBs of RBG. 'N' is the number of eligible
            %                        UEs and 'P' is the RBG size in RBs
            %       mcsRBG         - MCS for eligible UEs based on the CQI values of the RBs
            %                        of RBG. This is a N-by-2 matrix where 'N' is number of
            %                        eligible UEs. For each eligible UE it contains, MCS
            %                        index (first column) and efficiency (bits/symbol
            %                        considering both Modulation and Coding scheme)
            %       pastDataRate   - Served data rate. Vector of N elements containing
            %                        historical served data rate to eligible UEs. 'N' is
            %                        the number of eligible UEs
            %       bufferStatus   - Buffer-Status of UEs. Vector of N elements where 'N'
            %                        is the number of eligible UEs, containing pending
            %                        buffer status for UEs
            %       ttiDur         - TTI duration in ms
            %       UEs            - RNTI of all the UEs (even the non-eligible ones for
            %                        this RBG)
            %       lastSelectedUE - The RNTI of the UE which was assigned the last
            %                        scheduled RBG
            %
            %   SELECTEDUE The UE (among the eligible ones) which gets assigned
            %   this particular resource block group
            %
            %   MCSINDEX The suitable MCS index based on the channel conditions

            selectedUE = -1;
            maxPFWeightage = 0;
            mcsIndex = -1;
            linkDir = schedulerInput.LinkDir;
            for i = 1:length(schedulerInput.eligibleUEs)
                bufferStatus = schedulerInput.bufferStatus(i);
                pastDataRate = obj.UEsServedDataRate(schedulerInput.eligibleUEs(i), linkDir+1);
                if(bufferStatus > 0) % Check if UE has any data pending
                    
%                     schedulerInput.cqiRBG
%                     schedulerInput.mcsRBG
                    bitsPerSym = schedulerInput.mcsRBG(i, 2); % Accounting for both Modulation & Coding scheme
                    achievableDataRate = ((schedulerInput.RBGSize * bitsPerSym * 14 * 12)*1000)/ ...
                        (schedulerInput.ttiDur); % bits/sec
                    % Calculate UE weightage as per PF strategy
                    pfWeightage = achievableDataRate/pastDataRate;
                    if(pfWeightage > maxPFWeightage)
                        % Update the UE with maximum weightage
                        maxPFWeightage = pfWeightage;
                        selectedUE = schedulerInput.eligibleUEs(i);
                        mcsIndex = schedulerInput.mcsRBG(i, 1);
                    end
                end
            end

        end

    end
    methods(Access = protected)
            
            % Update 
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
                    % Update assignULResourceTTI
                    uplinkGrants = assignULResourceTTIConfigured2(obj, slotNum, firstULSym, numULSym);

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
        function downlinkGrants = scheduleDLResourcesSlot(obj, slotNum)
            %scheduleDLResourcesSlot Schedule DL resources of a slot
            % Downlink grants are returned as output to convey the way the
            % the downlink scheduler has distributed the resources to
            % different UEs. 'slotNum' is the slot number in the 10 ms
            % frame which is getting scheduled. The output 'downlinkGrants' is
            % a cell array where each cell-element represents a downlink
            % grant and has following fields:
            %
            % RNTI        Downlink grant is for this UE
            %
            % Type        Whether assignment is for new transmission ('newTx'),
            %             retransmission ('reTx')
            %
            % HARQId   Selected downlink UE HARQ process ID
            %
            % RBGAllocationBitmap  Frequency domain resource assignment. A
            %                      bitmap of resource block groups of the PUSCH
            %                      bandwidth. Value 1 indicates RBG is assigned
            %                      to the UE
            %
            % StartSymbol  Start symbol of time-domain resources. Assumed to be
            %              0 as time-domain assignment granularity is kept as
            %              full slot
            %
            % NumSymbols   Number of symbols allotted in time-domain
            %
            % SlotOffset   Slot-offset of PUSCH assignments for upcoming slot
            %              w.r.t the current slot
            %
            % MCS          Selected modulation and coding scheme for UE with
            %              respect to the resource assignment done
            %
            % NDI          New data indicator flag
            %
            % FeedbackSlotOffset Slot offset of PDSCH ACK/NACK from PDSCH transmission (i.e. k1)

            % Calculate offset of the slot to be scheduled, from the current
            % slot
            if slotNum >= obj.CurrSlot
                slotOffset = slotNum - obj.CurrSlot;
            else
                slotOffset = (obj.NumSlotsFrame + slotNum) - obj.CurrSlot;
            end

            % Get start DL symbol and number of DL symbols in the slot
            if obj.DuplexMode == 1 % TDD mode
                DLULPatternIndex = mod(obj.CurrDLULSlotIndex + slotOffset, obj.NumDLULPatternSlots);
                slotFormat = obj.DLULSlotFormat(DLULPatternIndex + 1, :);
                firstDLSym = find(slotFormat == obj.DLType, 1, 'first') - 1; % Location of first DL symbol in the slot
                lastDLSym = find(slotFormat == obj.DLType, 1, 'last') - 1; % Location of last DL symbol in the slot
                numDLSym = lastDLSym - firstDLSym + 1;
            else
                % For FDD, all symbols are DL symbols
                firstDLSym = 0;
                numDLSym = 14;
            end

            if obj.SchedulingType == 0  % Slot based scheduling
                % Assignments to span all the symbols in the slot
                downlinkGrants = assignDLResourceTTI(obj, slotNum, firstDLSym, numDLSym);
                % Update served data rate for the UEs as per the resource
                % assignments. This affects scheduling decisions for future
                % TTI
                updateUEServedDataRate(obj, obj.DLType, downlinkGrants);
            else %Symbol based scheduling
                if numDLSym < obj.TTIGranularity
                    downlinkGrants = [];
                    return; % Not enough symbols for minimum TTI granularity
                end
                numTTIs = floor(numDLSym / obj.TTIGranularity); % DL TTIs in the slot

                % DL grant array with maximum size to store grants. Maximum
                % grants possible in a slot is the product of number of
                % TTIs in slot and number of UEs
                downlinkGrants = cell((ceil(14/obj.TTIGranularity) * length(obj.UEs)), 1);
                numDLGrants = 0;

                % Schedule all DL TTIs in the slot one-by-one
                startSym = firstDLSym;
                for i = 1 : numTTIs
                    TTIDLGrants = assignDLResourceTTI(obj, slotNum, startSym, obj.TTIGranularity);
                    downlinkGrants(numDLGrants + 1 : numDLGrants + length(TTIDLGrants)) = TTIDLGrants(:);
                    numDLGrants = numDLGrants + length(TTIDLGrants);
                    startSym = startSym + obj.TTIGranularity;

                    % Update served data rate for the UEs as per the resource
                    % assignments. This affects scheduling decisions for future
                    % TTI
                    updateUEServedDataRate(obj, obj.DLType, TTIDLGrants);
                end
                downlinkGrants = downlinkGrants(1 : numDLGrants);
            end
        end
        
       % Change
       function eligibleUEs = getEligibleUEs(obj)
           eligibleUEs= [1,2,3,4,5,6,7,8,9,10];
           for n=1: height(obj.configurations)
           eligibleUEs(eligibleUEs == obj.configurations(n,1)) = [];
           end
        end

    end

    methods(Access = private)
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
        % Change 
            function ULGrantsTTI = assignULResourceTTIConfigured2(obj, slotNum, startSym, numSym)
            %assignULResourceTTI Perform the uplink scheduling of a set of contiguous UL symbols representing a TTI, of the specified slot
            % A UE getting retransmission opportunity in the TTI is not
            % eligible for getting resources for new transmission. An
            % uplink assignment can be non-contiguous, scattered over RBGs
            % of the PUSCH bandwidth
            
            grant= configuredULGrant(obj,slotNum);
            if ~isempty(grant)
            combinedRBGAllocationBitmap = combineRBG(obj,grant);
     
            end 
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

        % Change 
        function  releaseGrant(obj,RNTI,slotNum)
            for i=1:length(obj.ConfiguredGrants)
               if obj.ConfiguredGrants(i,1)==RNTI && obj.ConfiguredGrants(i,2)==slotNum
                   y(i)=false;
               else
                   y(i)=true;
               end
                
            end
            obj.ConfiguredGrants=obj.ConfiguredGrants(y,:);
        end
        
        % Change 
        function combinedRBG =combineRBG(obj,grant)
            if(~isempty(fieldnames(grant)))
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
            combinedRBG=zeros(1, obj.NumRBGsUL);
        end

        % Change
        function grant =configuredULGrant(obj,slot)
            
           grant = struct();            
           for n=1: height(obj.ConfiguredGrants)
            PRBs=  obj.ConfiguredGrants(n,2);
            UE= obj.ConfiguredGrants(find(obj.ConfiguredGrants(:,2)==PRBs),1);
              grant(n).RNTI = UE;  
            grant(n).Type = 'newTx';
            
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
                grant(n).RBGAllocationBitmap = [1 1 0 1 1 0 0 1 0 0 0 0 0];
            elseif (PRBs==7)
                grant(n).RBGAllocationBitmap = [0 1 1 1 1 1 1 1 0 0 0 0 0];
            elseif (PRBs==8)
                grant(n).RBGAllocationBitmap = [1 1 1 1 1 1 1 1 1 1 1 1 1];
            else
                grant(n).RBGAllocationBitmap = [0 0 0 0 0 0 0 0 0 0 0 0 0];
            end
            grant(n).StartSymbol = 0;
            grant(n).NumSymbols = 14;
            grant(n).SlotOffset = 1;
            grant(n).MCS = 1;
            grant(n).RV=0;
            grant(n).HARQId = 0; 
             if obj.prev==1 && obj.counter==3
                 obj.counter=1;
                 grant(n).NDI = 0; 
            	 obj.prev=0;
                 grant(n).Type = 'reTx';
                 
             else 
                 grant(n).NDI = 1; 
                  obj.prev=1;
                  obj.counter=obj.counter+1;
                  
             end 
            
            end 
       
 

    end
    end
end