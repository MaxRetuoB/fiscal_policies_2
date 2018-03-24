/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lars Ljunqsvist and Thomas Sargent model of fiscal policies in growth 
model from chapter 11  of "Recursive macroeconomic theory" 3rd ed.

Let's reproduce the figure 11.9.3 chap11, p398: "response to PERMANENT increase in 
g(governement spendingin g at t=10, and yield curve

   for this script we use two several equations:

                qt=((beta^(t))*(c(t))^(-gamma))/(c(0))^(-gamma)


                yit=-(1/i)*log(q(t+i)/q(t))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

Remark : the script model is available on the internet via "Practicing Dynare" 
(Barillas, Bhandari, Colacito, Kitao, Matthes, Sargent, Shin,2010), however
it seems that the script does not work with dynare 4.5.3. So it needs to be adapted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* the model:
Sargent and Ljungvist first model is characterised by an inelastic labor supply,
an exogenous stream of government expenditures and several kinds of distoring taxes.
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Declares the endogenous variables;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var c k;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//declares the exogenous variables:

//consumption tax (tauc), 
//capital tax (tauk), 
// government spending (g) 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	varexo tauc tauk g;

	parameters bet gam del alpha A;
	bet=.95; // discount factor
	gam=2; // CRRA parameter
	del=.2; // depreciation rate
	alpha=.33; // capital's share
	A=1; // total factor productivity
	/*------------------------------------------------------------
    model
    ------------------------------------------------------------*/
	model;
    // equation 11.3.8a:

	k=A*k(-1)^alpha+(1-del)*k(-1)-c-g;

    //equation 1.3.8e + 11.3.8.g
	c^(-gam) = bet*(c(+1)^(-gam))*((1+tauc)/(1+tauc(+1)))*((1-del)+
	(1-tauk(+1))*alpha*A*k(-1)^(alpha-1));
	end;
	
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Computation of the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    initval;
	k=1.5;
	c=0.6;
	g = .2;
	tauc = 0;
	tauk = 0;
	end;
	steady;

// The following values determine the new steady state
// after the shocks:
	
    endval;
	k=1.5;
	c=0.6;
	g =.4;
	tauc =0;
	tauk =0;
	end;
    steady;
	
// The following lines produce a g sequence with a once and for all jump in g shocks
// we use shocks to undo that for the first 9 periods and leave g at
// it?s initial value of 0.2
	
	
   

    shocks;
    var g;
    periods 1:9;
    values 0.2;
    end;

    simul(periods=1000);  

	 /*---------------------------------------------------------------
    graphs and plots for other endogenous variables
    ---------------------------------------------------------------*/
    
    //compute the initial steady state for consumption and capital for the plots
    c0=c(1);
    k0=k(1);
    g0=ex0_(3,1);

    // let N be the periods to plot
    N=40;

    /*-----------------------------------------------------------------
    the following equations compute the other endogenous variables use 
    in the plots below. These equations were taken from RMT 3 p 391
    -----------------------------------------------------------------*/

    // equation 11.6.8h

    rbig0=1/bet;

    rbig=c(2:1001).^(-gam)./(bet*c(3:1002).^(-gam));
    

    //equation 11.6.8c

    nuq0=alpha*A*k0^(alpha-1);

    nuq=alpha*A*k(1:1000).^(alpha-1);

    // equation 11.6.8d

    wq0=A*k0^alpha-k0*alpha*A*k0^(alpha-1);

    wq=A*k(1:1000).^alpha-k0*alpha*A*k(1:1000).^(alpha-1);

    // computing the price system: q and steady sttate price system
    // qs(t)=beta^t
    // q(2) is normalized to 1
    //qs(1) is normalized to 
    for n=1:102
    q(n)=(bet^n)*c(n).^(-gam)/((bet^2)*c(2)^(-gam));
    qs(n)=(bet^n)*c(1).^(-gam)/(bet*c(1)^(-gam));
    end

    // computing the short term interest rate
    rsmall=-log(q(3:102)./q(2:101));
    rsmall0=rbig0-1;

    //computing the term structure for n=1 to 40 for all periods
    for t=1:62
    for i=1:40
    y(i,t)=log(q(t+i)/q(t))/(-i);
    end
    end

    /*------------------------------------------------------------------
    Now we plot the yield curves and other responses 
    of the endogenous variables to the shock
    ------------------------------------------------------------------*/
    
    figure(1)
    x=0:N-1;
    

    //subplot for consumtpion 'c'
    subplot(2,3,1)
    plot(x,c0*ones(N,1),'--k',x,c(2:N+1),'k','Linewidth',1.5)
    title('consumption c','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)

   
    //subplot for consumption 'q'
     subplot(2,3,2)
    plot(x,qs(1:N),'--k',x,q(2:N+1),'k','Linewidth',1.5)
    title('consumption / investment price q','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)

    //subplot for rate of interest r
     subplot(2,3,3)
    plot(x,[rsmall0*ones(N,1)],'--k',x,rsmall(1:N),'k','Linewidth',1.5)
    //note the timing we lag capital to correct ofr syntax
    title('rate of interest $r$','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)

    //subplot for yield curves at t=0, t=10, t=60'
     subplot(2,3,4)
    plot(x,y(:,2),'k',x,y(:,11),'-.k',x,y(:,61),'--k','Linewidth',1.5)
    //note the timing we lag capital to correct ofr syntax
    title('yield curves at t=0, t=10, t=60','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)

    //subplot for the experiment proposed'
    subplot(2,3,5)
    plot([0:9],oo_.exo_simul(1:10,3),'k','Linewidth',1.5);
    hold on;
    plot([10:N-1],oo_.exo_simul(11:N,3),'k','Linewidth',1.5);
    hold on;
    plot(x,[g0*ones(N,1)],'--k','LineWidth',1.5)
    title('Government spending','interpreter','latex','Fontsize',12)
    axis([0 N -.1 .5])
    set(gca,'Fontsize',12)

    print -depsc fig_g_yields.eps
   
