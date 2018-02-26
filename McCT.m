function[Ntray,Rmin,L,V,Lbar,Vbar,B,D]=McCT(F,x1f,x1d,x1b,Rratio,q)


figure2=figure('Color',[1 1 1]);
axes2=axes('Parent',figure2,'FontSize',12);
box(axes2,'on')
hold(axes2,'all')  


  load('TempData');y_data=XY(:,2);x_data=XY(:,1);
     plot(x_data,y_data,'*')
    hold on
    plot([0 1],[0 1],'-k')
    hold on
    
    
    

     syms x f

    beta0=[-0.87 -0.1];
    modelfun=@(b,x)b(1)*exp(b(2)*x).*x.^2+(1-b(1)*exp(b(2)*x)).*x(:,1);
    mdl=fitnlm(XY(:,1),XY(:,2),modelfun,beta0);
    theta=mdl.Coefficients.Estimate;
    f=theta(1)*exp(theta(2)*x)*x^2+(1-theta(1)*exp(theta(2)*x))*x;
    warning off;
    %Will get 'cannot solve symbolically; returning a numeric approximation
    %instead'
    ezplot(f,[0 1])
    title([])
    hold on
    
    if q==1
        q=q+0.00001;
    end; %slight perturbation to avoid error associated with infinity slope
    
    %q line
    syms qline
    qline=q*x/(q-1)-x1f/(q-1);
    eqn=qline-f==0;
    solx=double(solve(eqn,x));
    solx=choice(solx);                          % function call takes result that is between
                              % 0 and 1
      ezplot(qline,[x1f solx]); 
    title([])
    hold on
    Rmin=x1d-(x1d-subs(qline,solx))*x1d/(x1d-solx);
    Rmin=x1d/Rmin-1;
    Rmin=double(Rmin);
    R=Rmin*Rratio;

    L_V=R/(R+1);
    
    DBvec=inv([1 1;x1d x1b])*[F;F*x1f];
    D=DBvec(1);B=DBvec(2);
    V=D/(1-L_V);
    L=V*L_V;
    Lbar=L+q*F;
    Vbar=V-(1-q)*F;
    
    syms enrichline
    enrichline=R*x/(R+1)+x1d/(R+1);
    ezplot(enrichline,[0 1]);
    title([])
    
    %Find intersection of q-line and enriching operation line
    eqn=qline-enrichline==0;
    solx=double(solve(eqn,x));
    x0=x1b;
    y0=x1b;
    x1=solx;
    y1=double(subs(qline,solx));
    syms stripline %Bottom line
    %stripline=y0+(y1-y0)*(x-x0)/(x1-x0);
    stripline=(Lbar/Vbar)*x-(Lbar/Vbar-1)*x1b;
    ezplot(stripline,[x0 1]);
    title([])
    
    %Now start adding in stages
    xact=x1d;
    actline=enrichline;
    iter=0;
    
    while xact>x1b
    yval2=double(subs(actline,xact));
    x2=double(solve(f==yval2,x));
    if length(x2)>1
    x2=choice(x2);
    end
    plot([xact x2],[yval2 yval2],'-k')
    hold on
    xact=x2;
    topop=double(subs(enrichline,xact));
    topop=topop(topop<1);
    botop=double(subs(stripline,xact));
    botop=botop(botop<1);
    yval1=min([topop botop]);
    if botop<topop
        actline=stripline;
    end
    plot([xact xact],[yval1 yval2],'-k')
    iter=iter+1;
    if iter>100
        break
    end
    end
    Ntray=iter-1;
    xlim([0 1])
    ylim([0 1])
    xlabel('x_a','FontSize',12,'FontWeight','bold')
ylabel('y_a','FontSize',12,'FontWeight','bold')
end
