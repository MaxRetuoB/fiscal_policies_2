/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lars Ljunqsvist and Thomas Sargent model of fiscal policies in growth 
model from chapter 11 of "Recursive macroeconomic theory" 3rd ed.

Know we make the assumption that the household as an elastic supply labor 

Let's reproduce the figure 11.12.1 chap11, p416: "Elastic labor suplly
response to PERMANENT increase in g(governement spending) at t=0
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

var c k n;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//declares the exogenous variables:

//consumption tax (tauc), 
//capital tax (tauk), 
// government spending (g) 
// income tax (taun)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	varexo tauc tauk taun g;

	parameters bet gam del alpha A B;
	bet=.95; // discount factor
	gam=1; // CRRA parameter
	del=.2; // depreciation rate
	alpha=.33; // capital's share
	A=1; // total factor productivity
    B=3; //coefficient on leisure

	/*------------------------------------------------------------
    model
    ------------------------------------------------------------*/
	
    model;
    //feasibility 
   
	k=A*k(-1)^alpha*n^(1-alpha)+(1-del)*k(-1)-c-g;

   
    //Euler equation
	c^(-gam) = bet*(c(+1)^(-gam))*((1+tauc)/(1+tauc(+1)))*((1-del)+
	(1-tauk(+1))*alpha*A*k^(alpha-1)*n(+1)^(1-alpha)+tauk(+1)*del);

    //consumption leisure choice
    B/c^(-gam)=((1-taun)/(1+tauc))*((1-alpha)*A*k(-1)^(alpha)*n^(-alpha));
	end;
	
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Computation of the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    initval;
	k=1.5;
	c=0.6;
	g = .2;
    n=1.02;
	tauc = 0;
    taun =0;
	tauk = 0;
	end;
	steady;

// The following values determine the new steady state
// after the shocks:
	
    endval;
	k=1.5;
	c=0.6;
	g =.4;
    n=1.02;
	tauc =0;
    taun =0;
	tauk =0;
	end;
    steady;
	
// The following lines produce a g sequence with a once and for all jump in g shocks
// we use shocks to undo that for the first 9 periods and leave g at
// it?s initial value of 0.2
	
	
   /*

    shocks;
    var g;
    periods 1:9;
    values 0.2;
    end;
*/
    simul(periods=1000);  

	 /*---------------------------------------------------------------
    graphs and plots for other endogenous variables
    ---------------------------------------------------------------*/
    
    //compute the initial steady state for consumption and capital for the plots
    
     c0=c(1);
    k0=k(1);
    n0=n(1)
    g0=ex0_(4,1);
   
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
    
    wq0=(1-alpha)*A*(k0^alpha)*(n0^(-alpha));

    wq=(1-alpha)*A*(k(1:100).^alpha).*(n(2:101).^(-alpha));
    y=A*(k(1:100).^alpha).*(n(2:101).^(1-alpha));
    gr=y(2:100)-y(1:99);
    grf=gr(1:99)/y(1:99);
    /*------------------------------------------------------------------
    Now we plot the responses of the endogenous variables to the shock
    ------------------------------------------------------------------*/
    x=0:N-1;
    figure(2)
    plot(x,gr(1:N),'k','Linewidth',1.5)
    title('Economic growth','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)
    x=0:N-1;

    x=0:N-1;
    figure(1)
     
    //subplot for capital 'k'
    
    subplot(2,3,1)
    plot(x,[k0*ones(N,1)],'--k',x,k(1:N),'k','Linewidth',1.5)
    
    //note the timing we lag capital to correct ofr syntax
    
    title('capital','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)
    
   
    //subplot for consumption
    
     subplot(2,3,2)
    plot(x,[c0*ones(N,1)],'--k',x,c(2:N+1),'k','Linewidth',1.5)
    title('consumption','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)

    //subplot for leisure
    
    subplot(2,3,3)
    plot(x,[n0*ones(N,1)],'--k',x,n(2:N+1),'k','Linewidth',1.5)
    title('Labor','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)
    
    //subplot for cost of capital 'Rbarr'
    
     subplot(2,3,4)
    plot(x,[rbig0*ones(N,1)],'--k',x,rbig(1:N),'k','Linewidth',1.5)
    title('Cost of capital $\overline{R}$','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)
    
    //subplot for rental rate 'wage'
    
     subplot(2,3,5)
    plot(x,[wq0*ones(N,1)],'--k',x,wq(1:N),'k','Linewidth',1.5)
    title('wage','interpreter','latex','Fontsize',12)
    set(gca,'Fontsize',12)
    
    //subplot for the experiment proposed'
    
    subplot(2,3,6)
    plot([0:40],oo_.exo_simul(2:42,4),'k','Linewidth',1.5)
    hold on;
    
    plot(x,[g0*ones(N,1)],'--k','LineWidth',1.5)
    title('Government spending','interpreter','latex','Fontsize',12)
    axis([0 N -.1 .5])
    set(gca,'Fontsize',12)

    print -depsc fig_g_elas.eps
   
