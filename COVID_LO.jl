using PyPlot

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

Infetti = [15,55,112,173,240,305,403,531,615,984,1254];
Morti = [0,1,2,6,8,9,14,17,23,24,38];
Guariti = [0,0,0,0,0,0,42,46,50,83,149];

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
	X,Y = EulerMethod(dy,[S0,I0,R0],[0.0,0.0,0.0],[0.01,0.01,0.01],160*(1/0.01));
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
	I0 = 1.0;
	R0 = 0.0;
	S0 = N-I0-R0;
	dy(y) = [ -β*y[1]*(y[2]/N),
		  β*y[1]*(y[2]/N) - γ*y[2],
		  γ*y[2]
		];
	X,Y = EulerMethod(dy,[S0,I0,R0],[0.0,0.0,0.0],[0.01,0.01,0.01],d*(1/0.01));
	#title(L"Closed SIR Model $\beta=0.2$, $\gamma=0.1$");
	#plot(X,Y);
	#semilogygend(["Susceptible","Infected","Recoverd and Immune"]);
	return X,Y;
end
function LSE(Y,D1,D2)
	S = 0
	for i in 1:length(D1)
		S = S + ((10^-5)*(Y[i][2]-D1[i]))^2+((10^-5)*(Y[i][3]-D2[i]))^2;
		#println(S);
	end
	return S;

end
function GradientDescent(itmax,subit,h)
	β = 0.8; γ = 0.10; w=10;
	E = [1000000000.0];
	L = 1000000000.0;
	for it in 1:itmax
		for i in 1:subit
			X,Y = SIR_Closed(β,γ,11);
			L = LSE(Y,Infetti,Morti.+Guariti);
			#println("iit: ", i, "L: ", L);
			X,Y = SIR_Closed(β-h,γ,11);
			L1 = LSE(Y,Infetti,Morti.+Guariti);
			X,Y = SIR_Closed(β+h,γ,11);
			L2 = LSE(Y,Infetti,Morti.+Guariti);
			∂β = (L2-L1)/(2*h);  
			X,Y = SIR_Closed(β,γ-h,11);
			L1 = LSE(Y,Infetti,Morti.+Guariti);
			X,Y = SIR_Closed(β,γ+h,11);
			L2 = LSE(Y,Infetti,Morti.+Guariti);
			∂γ = (L2-L1)/(2*h);
			β = β-w*∂β;
			γ = γ-w*∂γ;
			if E[end] < L
				w = 0.5*w;
			else
				break
			end
		end
		println("Iteration: ",it," step:", h," Energy:", L);
		E = push!(E,L);

	end
	plot(E);
	return β,γ;
	
end
function SIR_Closed_Fitted()

	β,γ = GradientDescent(1000,100,0.001);
	println("β=",β," γ=",γ);
	figure()
	X,Y = SIR_Closed(0.8,1.0/10,11);
	title("Closed SIR Model");
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	scatter(1:11,Infetti)
	scatter(1:11,Morti.+Guariti)
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Infected","Removed","True Infected","True Removed"]);

	figure()
	X,Y = SIR_Closed(β,γ,11);
	title("Closed SIR Model");
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	scatter(1:11,Infetti)
	scatter(1:11,Morti.+Guariti)
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Infected","Removed","True Infected","True Removed"]);

	figure()
	X,Y = SIR_Closed(β,γ,360);
	title("Closed SIR Model");
	plot([x[1] for x in X],[y[1] for y in Y]);
	plot([x[2] for x in X],[y[2] for y in Y]);
	plot([x[3] for x in X],[y[3] for y in Y]);
	scatter(1:11,Infetti)
	scatter(1:11,Morti.+Guariti)
	println("Least square error: ",LSE(Y,Infetti, Morti.+Guariti));
	legend(["Susceptible","Infected","Removed","True Infected","True Removed"]);
	show();
end
#Euler1DEx();
#LotkaVolterraEx();
#SIR_Closed();
SIR_Closed_Fitted();
