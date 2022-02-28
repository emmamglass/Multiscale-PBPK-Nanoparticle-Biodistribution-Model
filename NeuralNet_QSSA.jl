using DifferentialEquations
using Plots
using Flux
using Statistics

tmax = 1000.0
names = ["lungs", "heart", "kidney", "liver", "spleen",  "Gut tissue", "Other Tissue"]  

#Model Parameters
A_i=[.0019956689, .0009621096822, .0033612343, .0045212562, .0004750356, .003333333, .003333333].*1000 #mouse dm^2 
l_NC_b  = (100*10^(-9)) # m
Q_i=[74.83758415, 138.673326, 763.9168803, 284.6716878, 224.48828292, 333, 1000]./(60*10^6) #L/sec MOUSE
V_i_BL=[43.65525742 ,20.04395171 ,1018.55584 ,127.7673516 ,2.226207538 ,90 ,700]./10^6
V_i_T=[33.26114851 ,16.03516137 ,56.02057122 ,75.3542703 ,7.91726 ,20 ,1000]./10^6 #L MOUSE

V_vein = 466.9/10^6 # mouse L
Q_hep = sum(Q_i[4:6]) #flow out of the liver
V_hep = sum(V_i_BL[4:6]) #total volume of blood passed thru liver
Q_vein = Q_i[3]+Q_hep+Q_i[7] #L/sec Mouse

#Rate Constants
#Kon and Koff
slope1 = (22-4)/(500-50)
slope2 = 3.9
slope = slope1*slope2
a = [4, 15, 50, 79, 100]
r = a/2
K_EC0=[1.23e42, 10.6565662, 41793.3924, 2.28e12, 7.74e21, 3.33e15, 3.33e15]
K_EC0 = log.(K_EC0)
K_EC = zeros(5, 7)
for i =1:7
    K_EC[1:4,i] = slope*r[1:4] .- slope*r[end] .+ K_EC0[i]
    K_EC[end, i] = K_EC0[i]
end
K_EC[:, 2] .= K_EC0[2]

Kb = 1.36e-23
T = 310
mu = 0.0037
l = 1e-7
r = r*1e-9
D=Kb*T./(6*pi*mu*r)
k_i_on = D./l^2
K_i_on = zeros(5, 7)
K_i_off = zeros(5, 7)
for i = 1:5
  K_i_on[i, :] .= k_i_on[i]
  K_i_off[i, :] = K_i_on[i]./K_EC[i, :]  
end

K_i_on = K_i_on[1, :]
K_i_off = K_i_off[1, :]
K_deg_vein = 1/6000

#Other Parameters
p = [2.0e-6 1.0 1.0; 2e-6 0.1 0.02; 0.0004 0.0033 0.02; 5.0e-6 0.005 0.02; 2e-5 0.08 0.014; 3.3e-5 0.01 0.02; 3.3e-5 0.01 0.02]
#Rate of Degredation of NP in each organ
K_i_deg = p[:, 1]
#Rate of uptake of NP in each organ
K_i_up = p[:, 2]
#Rate of non specific uptake of NP in each organ
K_NS= p[:, 3]

#ODE matrix for linear system (Ay = b)
A = [-(Q_i[1]/V_i_BL[1]+K_i_on[1]+K_NS[1]+K_i_deg[1]) Q_i[1]/V_i_BL[2] 0 0 0 0 0 K_i_off[1] 0 0 0 0 0 0;
Q_i[1]/V_i_BL[1] -((Q_i[1]+Q_vein)/V_i_BL[2] + K_i_on[2]+K_NS[2]+K_i_deg[2]) 0 0 0 0 0 0 K_i_off[2] 0 0 0 0 0;
0 0 -(Q_i[3]/V_i_BL[3] + K_i_on[3]+K_NS[3]+K_i_deg[3]) 0 0 0 0 0 0 K_i_off[3] 0 0 0 0;
0 0 0 -(Q_hep/V_i_BL[4]+K_i_on[4]+K_NS[4]+K_i_deg[4]) Q_i[5]/V_i_BL[5] Q_i[6]/V_i_BL[6] 0 0 0 0 K_i_off[4] 0 0 0;
0 0 0 0 -(Q_i[5]/V_i_BL[5]+K_i_on[5]+K_NS[5]+K_i_deg[5]) 0 0 0 0 0 0 K_i_off[5] 0 0;
0 0 0 0 0 -(Q_i[6]/V_i_BL[6] + K_i_on[6] + K_i_deg[6] + K_NS[6]) 0 0 0 0 0 0 K_i_off[6] 0;
0 0 0 0 0 0 -(Q_i[7]/V_i_BL[7] + K_i_on[7] + K_NS[7] + K_i_deg[7]) 0 0 0 0 0 0 K_i_off[7];
K_i_on[1] 0 0 0 0 0 0 -(K_i_up[1] + K_i_off[1] + K_i_deg[1]) 0 0 0 0 0 0;
0 K_i_on[2] 0 0 0 0 0 0 -(K_i_up[2] + K_i_off[2] + K_i_deg[2]) 0 0 0 0 0;
0 0 K_i_on[3] 0 0 0 0 0 0 -(K_i_up[3] + K_i_off[3] + K_i_deg[3]) 0 0 0 0;
0 0 0 K_i_on[4] 0 0 0 0 0 0 -(K_i_up[4] + K_i_off[4] + K_i_deg[4]) 0 0 0;
0 0 0 0 K_i_on[5] 0 0 0 0 0 0 -(K_i_up[5] + K_i_off[5] + K_i_deg[5]) 0 0;
0 0 0 0 0 K_i_on[6] 0 0 0 0 0 0 -(K_i_up[6] + K_i_off[6] + K_i_deg[6]) 0;
0 0 0 0 0 0 K_i_on[7] 0 0 0 0 0 0 -(K_i_up[7] + K_i_off[7] + K_i_deg[7])]

#ODE function with quasi-steady-state component
function func_qss!(dydt, y, p, t)
    b_vec = [0; -Q_vein*y[1]/V_vein; -Q_i[3]*y[2]/V_vein;-Q_i[4]*y[2]/V_vein;
            -Q_i[5]*y[2]/V_vein;-Q_i[6]*y[2]/V_vein;
            -Q_i[7]*y[2]/V_vein;
             0;0;0;0;0;0;0]
  
    #Steady-State component
    y_ss = A\b_vec
  
    #Vein
    dydt[1] = (Q_i[3]*y_ss[3]/V_i_BL[3] + Q_hep*y_ss[4]/V_i_BL[4] + Q_i[7]*y_ss[7]/V_i_BL[7] - Q_vein*y[1]/V_vein - y[1]*K_deg_vein) 
    # Artery 
    dydt[2] = (Q_vein*y_ss[2]/V_i_BL[2]-(Q_i[3]*y[2] + Q_i[4]*y[2]+ Q_i[5]*y[2] + Q_i[6]*y[2] + Q_i[7]*y[2])/V_vein - y[2]*K_deg_vein) 
  
    # ---------------------- organ tissue ------------------------ #
    dydt[3:9] = (K_i_up.*y_ss[8:14] + K_NS.*y_ss[1:7]-K_i_deg[1:end].*y[3:9])
end

#Solving ODE with non-stiff solver
y0=[1.0, 0, 0, 0, 0, 0, 0, 0, 0]
tspan = (0, 1e3)
prob = ODEProblem(func_qss!, y0, tspan)
sol_qss = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10)

#Solving ODEs with QSSA with a Neural Network
function ODEfunc!(y)

    b_vec = [zeros(1, bat_sz); -Q_vein*y[1, :]'/V_vein; -Q_i[3]*y[2, :]'/V_vein;
            -Q_i[4]*y[2, :]'/V_vein;-Q_i[5]*y[2, :]'/V_vein;
            -Q_i[6]*y[2, :]'/V_vein;-Q_i[7]*y[2, :]'/V_vein;
            zeros(1, bat_sz);
            zeros(1, bat_sz);
            zeros(1, bat_sz);
            zeros(1, bat_sz);
            zeros(1, bat_sz);
            zeros(1, bat_sz);
            zeros(1, bat_sz)]
    
    #Steady State Component
    y_ss = A\b_vec
    a = (Q_i[3]*y_ss[3, :]/V_i_BL[3] .+ Q_hep*y_ss[4, :]/V_i_BL[4] .+ Q_i[7]*y_ss[7, :]/V_i_BL[7] 
        .- Q_vein*y[1, :]/V_vein .- y[1, :]*K_deg_vein)
    # Artery 
    b = (Q_vein*y_ss[2, :]/V_i_BL[2] .- (Q_i[3]*y[2, :] .+ Q_i[4]*y[2, :] .+ Q_i[5]*y[2, :] 
            .+ Q_i[6]*y[2, :] .+ Q_i[7]*y[2, :])/V_vein .- y[2, :]*K_deg_vein) 

    # ---------------------- organ tissue ------------------------ #
    c = (K_i_up.*y_ss[8:14, :] .+ K_NS.*y_ss[1:7, :] .- K_i_deg[1:end].*y[3:9, :])
    return [a'; b'; c]
              
end

function NeuralNetwork()
    return Chain(Dense(1, 32, tanh),
    Dense(32, 32, tanh),
    Dense(32, 32, tanh),
    Dense(32, 9, sigmoid) )
end

t = collect(0.0 : 0.0005 : 1.0)
bat_sz = length(t)
t = reshape(t, 1, :)
data = Flux.Data.DataLoader(t, batchsize=size(t)[2],shuffle=true)

m = NeuralNetwork()
# g(t) = t.*m(t) .+ y0
opt = Flux.ADAM(0.001)

ϵ = sqrt(eps(Float32))
y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#Derivative of output from neural network
derivative_fun(t) = (-m(t.+2ϵ)+8m(t.+ϵ)-8m(t.-ϵ)+m(t.-2ϵ))/(12ϵ)/tmax

#Loss function
loss(t) = Flux.Losses.mse(derivative_fun(t)  , ODEfunc!(m(t)); agg = mean) + Flux.Losses.mse(m([0.0]), y0; agg = mean)

#Neural Network Model Parameters
ps = Flux.params(m)
epochs = 30000
iter = 0

cb1 = function () #callback function to observe training
    global iter += 1
    if iter % 1000 == 0
        print(loss(t), "\n")
        flush(stdout)
    end
end

#Training Neural Network
for i in 1:epochs
    Flux.train!(loss, ps, data, opt, cb = cb1)
end
#Forward Pass through neural network
sol = m(t)

t = collect(0:0.5:1000)
y_ub = sol_qss(t)
y_val_ub = zeros(9 , size(y_ub.u)[1])
for i=1:size(y_ub.u)[1]
    for j=1:9
        y_val_ub[j, i] = y_ub.u[i][j]
    end
end
organs = ["NP (Lungs)", "NP (Heart)", "NP (Kidney)", "NP (Liver)", "NP (Spleen)", "NP (Gut)", "NP (Others)"]


#Neural Network are Tsit5 solution Plot together
plot_array = Any[]
gr()
push!(plot_array, plot([0.0001.+t/3600, 0.0001.+t/3600], [sol[1:2, :]', y_val_ub[1:2, :]'], xaxis=:log, label = ["Neural Net(Vein)" "Neural Net(Artery)" "ODE Solver(Vein)" "ODE Solver(Artery)"],ylabel="NP (Vein/Artery)",size=(800,800)))
for i=1:7
    push!(plot_array, plot([0.0001.+t/3600, 0.0001.+t/3600], [sol[i+2, :], y_val_ub[i+2, :]],xaxis=:log, ylabel=organs[i],label = ["Neural Net" "ODE solver"], legend=:bottomright,size=(800,800)))
end
display(plot(plot_array..., xlabel="Time (Hours)",layout=(4,2), legendfont=font(8)))
