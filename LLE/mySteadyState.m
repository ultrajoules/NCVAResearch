function FS = mySteadyState(uS,Ein,delta)
        FS = Ein./(1+1i*(delta - uS.*conj(uS)))-uS;
end