function [Tk,Fk,Ak,CFk,Xk]=label_vector(s,h)

Tk=zeros(size(s(:,1)));
Fk=zeros(size(s(:,1)));
Ak=zeros(size(s(:,1)));
CFk=zeros(size(s(:,1)));
Xk=zeros(size(s(:,1)));
count=0;

trial_start=1;
fixation_cross=786;
cue=[771,773,783];
continuous_feedback=781;
hit=897;
miss=898;

for idx=1:length(h.EVENT.TYP)
    
    if h.EVENT.TYP(idx)==trial_start
        count=count+1;
        Tk(h.EVENT.POS(idx):h.EVENT.POS(idx)+h.EVENT.DUR(idx))=count;
        i=1;
        while idx+i<=length(h.EVENT.TYP) & h.EVENT.TYP(idx+i)~=0 
            Tk(h.EVENT.POS(idx+i):h.EVENT.POS(idx+i)+h.EVENT.DUR(idx+i))=count;
            i=i+1;
        end
    end
    
    if h.EVENT.TYP(idx)==fixation_cross
        Fk(h.EVENT.POS(idx):h.EVENT.POS(idx)+h.EVENT.DUR(idx)-1)=fixation_cross;
    end

    if find(cue==h.EVENT.TYP(idx))
        Ak(h.EVENT.POS(idx):h.EVENT.POS(idx)+h.EVENT.DUR(idx)-1)=h.EVENT.TYP(idx);
    end

    if h.EVENT.TYP(idx)==continuous_feedback
        CFk(h.EVENT.POS(idx):h.EVENT.POS(idx)+h.EVENT.DUR(idx)-1)=continuous_feedback;
    end

    if h.EVENT.TYP(idx)==hit
        Xk(h.EVENT.POS(idx):end)=hit;
    end
    if h.EVENT.TYP(idx)==miss
        Xk(h.EVENT.POS(idx):end)=miss;
    end

end