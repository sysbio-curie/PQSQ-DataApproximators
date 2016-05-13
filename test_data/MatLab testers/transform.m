function res = transform( str )
%transform transforms the results of tests to readable form

    for k=1:length(str);
        st.name = str(k).name;
        st.errMean = mean(str(k).error);
        st.errStdev = std(str(k).error);
        st.timeMean = mean(str(k).time);
        st.timeStdev = std(str(k).time);
        if k==1
            res = st;
        else
            res = [res; st];
        end
    end
end

