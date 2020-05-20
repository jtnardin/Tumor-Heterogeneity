function HandlePropertyEvents(src,evnt)
    switch src.Name 
        case 'Kbar'
            % Kbar is changed, re-randomize centers
            evnt.AffectedObject.Randomize();
        case 'cov'
            % cov is changed, make sure it's the right size
            if(length(evnt.AffectedObject.cov)==1)
                % single lump width
            evnt.AffectedObject.cov = evnt.AffectedObject.cov*eye(evnt.AffectedObject.dim);  
            end
            if(all(size(evnt.AffectedObject.cov)==[evnt.AffectedObject.dim,1]))  
                % Anisotropic-but-diagonal
            evnt.AffectedObject.cov = diag(evnt.AffectedObject.cov);  
            end
            % Need to re-randomize centers because the edge padding is
            % now different
            if(evnt.AffectedObject.showwarnings)
                warning('The lump radius has been changed; edge effects may occur! Recommend re-randomizing.');
            end
            %evnt.AffectedObject.randomize();
        case 'b0'
            % b0 is changed, make sure it's the right size
            if(length(evnt.AffectedObject.b0)~=evnt.AffectedObject.K)
                evnt.AffectedObject.b0 = evnt.AffectedObject.b0*ones(evnt.AffectedObject.K,1);
            end
    end
end