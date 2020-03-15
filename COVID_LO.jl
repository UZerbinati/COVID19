using PyPlot
using DataFrames
using LsqFit

hE = 0.01;
function EulerMethod(dy,I,D,h,it)
	"""
		Implementation of the farward Euler method, for 1D system and ND systems.
	"""
	Y = [I]; #Defining initial condition.
	X = [D]; #Adding starting points.
	dx = dy(I).*h; #Computing first step
	for i in 1:it
		X = push!(X, X[end]+h);
		Y = push!(Y,Y[end]+dx); #Adding step
		dx = dy(Y[end]).*h; #Computing Next Step
	end
	return X,Y;
end

function Euler1DEx()
	"""
		Here we have the example of using the (farward) Euler method to compute the solution of dy(y) = y, with y=y(t),  wich has as solution exp(t).
	"""
	dy(y) = y;
	X,Y=EulerMethod(dy,1.0,0.0,0.01,100);
	plot(X,Y);
	plot(X,exp.(X));
	title("1D Euler Example");
	show();
end

function LotkaVolterraEx()
	"""
		Here we have the example of using the (farward) Euler method to compute the solution of the Lotka-Volterra differential equations.
	"""
	A = 1.5; B=1; C=3; D=1;
	dy(y) = [A*y[1] - B*y[1]*y[2], -C*y[2]+D*y[1]*y[2]];
	X,Y=EulerMethod(dy,[10.0, 5.0],[0.0,0.0],[0.01,0.01],1000);
	plot(X,Y);
	title("nD Euler Example");
	show();
end

Infetti = [15,55,112,173,240,305,403,531,615,984,1254,1520];
Morti = [0,1,2,6,8,9,14,17,23,24,38,55];
Guariti = [0,0,0,0,0,0,42,46,50,83,149,160];
LOBeta = [0.76,0.76,0.566,0.178,0.150,0.142];
function beta_SQ_fit()
	ω = log.(LOBeta);
	x = 1:4:23;
	#println(length(ω))
	#println(length(1:4:23))
	b = (sum(x.*ω)-(1/6)*sum(x)*sum(ω))/(sum(x.*x)-(1/6)*(sum(x)^2));
	a = (1/6)*sum(ω) - (1/6)*b*sum(x);
	d = 1:0.1:23;
	g(x) = exp(a)*exp(b*x);
	figure()
	scatter(x,LOBeta)
	plot(d,g.(d));
	title("β In Lodi Province")
	ylabel("β");
	xlabel("Day From The 21st Of February");
	legend(["β Forecast", "True β 4 Day Mean"]);
	return exp(a), b;
	
end
function SIR_Closed_Ex1()
	"""
		S'(t) = -βS(t)(I(t) / N)
		I'(t) = βS(t) ( I(t) / N ) - γ I(t)
		R'(t) = γ I(t)
	"""
	β=0.2;
	γ=1.0/10;
	N=1000.0;
	I0 = 1.0;
	R0 = 0.0;
	S0 = N-I0-R0;
	dy(y) = [ -β*y[1]*(y[2]/N),
		  β*y[1]*(y[2]/N) - γ*y[2],
		  γ*y[2]
		];
	X,Y = EulerMethod(dy,[S0,I0,R0],[0.0,0.0,0.0],[hE,hE,hE],160*(1/hE));
	title(L"Closed SIR Model $\beta=0.2$, $\gamma=0.1$");
	plot(X,Y);
	legend(["Susceptible","Infected","Recoverd and Immune"]);
	show();
end

function SIR_Closed(β,γ,d)
	"""
		S'(t) = -βS(t)(I(t) / N)
		I'(t) = βS(t) ( I(t) / N ) - γ I(t)
		R'(t) = γ I(t)
	"""
	N=10078012.0;
	I0 = 1;
	R0 = (Morti.+Guariti)[1];
	S0 = N-I0-R0;
	dy(y) = [ -β*y[1]*(y[2]/N),
		  β*y[1]*(y[2]/N) - γ*y[2],
		  γ*y[2]
		];
	X,Y = EulerMethod(dy,[S0,I0,R0],[0.0,0.0,0.0],[hE,hE,hE],d*(1/hE));
	#plot(X,Y);
	#legend(["Susceptible","Infected","Recoverd and Immune"]);
	return X,Y;
end
function SIR_Qu(β,γ,d,A,b,p)
	"""
		S'(t) = -βS(t)(I(t) / N)
		I'(t) = βS(t) ( I(t) / N ) - γ I(t)
		R'(t) = γ I(t)
	"""
	N=10078012.0;
	I0 = 1.0;
	R0 = 0.0;
	S0 = N-I0-R0;
	dy(y) = [ -(p*y[4]+(1.0-p))*β*y[1]*(y[2]/N),
		(p*y[4]+(1.0-p))*β*y[1]*(y[2]/N) - γ*y[2],
		  γ*y[2],
		  b*y[4]
		];
	X,Y = EulerMethod(dy,[S0,I0,R0,A],[0.0,0.0,0.0,0.0],[hE,hE,hE,hE],d*(1/hE));
	#plot(X,Y);
	#legend(["Susceptible","Infected","Recoverd and Immune"]);
	return X,Y;
end
function SIR_Test(β,γ)
	X,Y = SIR_Closed(β,γ,12);
	return LSE(Y,Infetti,Morti.+Guariti);
end
function SIR_Test_LSQ(β,γ)
	X,Y = SIR_Closed(β,γ,12);
	R = [];
	for i in 1:12
		push!(R,Y[floor(Int,(1/hE)*i)][2]);
	end
	return R;
end
function LSE(Y,D1,D2)
	S = 0
	for i in 2:length(D1)
		#println("Infetti: ", Y[floor(Int,(1/hE)*i)][2],"\tInfetti Reali: ", D1[i], "\tRimossi: ", Y[floor(Int,(1/hE)*i)][3], "\tRimossi Reali:", D2[i]);
		S = S +(1/D1[i])*abs((Y[floor(Int,(1/hE)*i)][2]-D1[i]))+(1/D2[i])*abs((Y[floor(Int,(1/hE)*i)][3]-D2[i]));
		#println(S);
	end
	return S;

end
function GD(β0,γ0,itmax,h)
	β=β0;
	w =10^(-5);
	E =  [];
	for it in 1:itmax			
		X,Y = SIR_Closed(β,γ0,12);
		L = LSE(Y,Infetti,Morti.+Guariti);
		println("Iteration  ",it,",\tEnergy: ", L,"\tBeta: ",β);
		E = push!(E,L);
		X,Y = SIR_Closed(β+h,γ0,12);
		L2 = LSE(Y,Infetti,Morti.+Guariti);
		X,Y = SIR_Closed(β-h,γ0,12);
		L1 = LSE(Y,Infetti,Morti.+Guariti);
		∂β = (L2-L1)/(2*h);
		β = β-w*∂β;
		X,Y = SIR_Closed(β,γ0,12);
		L = LSE(Y,Infetti,Morti.+Guariti);
		if E[end] <= L
			println("E: ", E[end]," L: ",L," flag: ", E[end]<=L);
			break;
		end

	end
	return β;
end
function GD2(β0,γ0,itmax,h)
	β=β0;
	γ=γ0;
	wβ =10^(-4);
	wγ = 10^(-4);
	E =  0;
	for it in 1:itmax			
		X,Y = SIR_Closed(β,γ,12);
		L = LSE(Y,Infetti,Morti.+Guariti);
		#println("Iteration  ",it,",\tEnergy: ", L,"\tBeta: ",β,"\t Gamma:",γ);
		E = L;
		X,Y = SIR_Closed(β+h,γ0,12);
		L2 = LSE(Y,Infetti,Morti.+Guariti);
		X,Y = SIR_Closed(β-h,γ0,12);
		L1 = LSE(Y,Infetti,Morti.+Guariti);
		∂β = (L2-L1)/(2*h);
		X,Y = SIR_Closed(β,γ+h,12);
		L2 = LSE(Y,Infetti,Morti.+Guariti);
		X,Y = SIR_Closed(β,γ-h,12);
		L1 = LSE(Y,Infetti,Morti.+Guariti);
		∂γ = (L2-L1)/(2*h);
		
		β = β-wβ*∂β;
		γ = γ-wγ*∂γ;
		X,Y = SIR_Closed(β,γ,12);
		L = LSE(Y,Infetti,Morti.+Guariti);
		if E <= L
			println("E: ", E," L: ",L," flag: ", E[end]<=L);
			break;
		end

	end
	return β,γ;
end
function SIR_Closed_Fitted(β,γ)
	
	figure()

	println("β=",β," γ=",γ);
	X,Y = SIR_Closed(β,γ,12);
	title(string("SIR Model 12 days β=",round(β,digits=3),", γ=",round(γ,digits=3)));
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	scatter(1:12,Infetti)
	scatter(1:12,Morti.+Guariti)
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Infected","Removed","True Infected","True Removed"]);

	figure()
	X,Y = SIR_Closed(β,γ,180);
	title(string("SIR Model 180 days β=",round(β,digits=3),", γ=",round(γ,digits=3)));
	plot([x[1] for x in X],[y[1] for y in Y]);
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Susceptible","Infected","Removed"]);

	println("R0: ", β/γ);
	show();
end

function GenerateTable()
	df = DataFrame(β = Float64[], γ = Float64[], E = Float64[],R0 = Float64[]);
	for β0 in 0.1:0.1:2
		println("Iteration Start, for β=",β0);
		for γ0 in 1/23:0.01:1
			β,γ=GD2(β0,γ0,10000,0.01);
			#β = GD(1.94,1/1.64,10000,0.01);
			L = SIR_Test(β,γ);
			push!(df,(β,γ,L,β/γ));
		end
		println("Iteration Done, for β=",β0);
	end
	return df;
end

#Euler1DEx();
#LotkaVolterraEx();
#SIR_Closed();
#SIR_Closed_Fitted(0.8,1/14);

function SIR_CG_Ex()
	
	figure()
	X,Y = SIR_Closed(0.76,1/6.8,12);
	title(string("SIR Model 12 days β=",0.76,", γ=",round(1/6.8,digits=3)))
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	scatter(1:12,Infetti)
	scatter(1:12,Morti.+Guariti)
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Infected","Removed","True Infected","True Removed"]);
	ylabel("Number Of People");
	xlabel("Day From The 21st Of February");

	figure()
	X,Y = SIR_Closed(0.76,1/6.8,120);
	title(string("SIR Model 120 days β=",0.76,", γ=",round(1/6.8,digits=3)))
	plot([x[1] for x in X],[y[1] for y in Y]);
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Susceptible","Infected","Removed","True Infected","True Removed"]);
	ylabel("Number Of People");
	xlabel("Day From The 21st Of February");
		


	show()
end
function SIR_Qu_Ex(P)
	A,b = beta_SQ_fit();
	println([A,b]);
	

	figure()
	for p=0.1:0.1:P
		X,Y = SIR_Qu(0.76,1/6.8,240,A,b,p);
		plot([x[2] for x in X],[y[2] for y in Y]);
	end
	title(string("SIR Model 240 days β=",0.76,", γ=",round(1/6.8,digits=3)))
	legend([string("Infected p=",p) for p in 0.1:0.1:P]);
	ylabel("Number Of People");
	xlabel("Day From The 21st Of February");
	
	show()
end
