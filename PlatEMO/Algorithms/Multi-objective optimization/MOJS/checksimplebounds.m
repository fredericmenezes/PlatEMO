function POS=checksimplebounds(POS,Lb,Ub)
    for i=1:size(POS,1)
        ns_tmp=POS(i,:);
        I=ns_tmp<Lb;
        while sum(I)~=0
            ns_tmp(I)=Ub(I)+(ns_tmp(I)-Lb(I));
            I=ns_tmp<Lb;
        end
        J=ns_tmp>Ub;
        while sum(J)~=0
            ns_tmp(J)=Lb(J)+(ns_tmp(J)-Ub(J));
            J=ns_tmp>Ub;
        end
        POS(i,:)=ns_tmp;
    end
end