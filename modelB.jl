using DifferentialEquations, LinearAlgebra
using Plots

names_new = ["vein", "artery", "Lung Vas", "Heart Vas", "Kidney Vas", "Liver Vas", "Spleen Vas", "Gut Vas" ,
        "Others Vas", "lungs et", "Heart et", "Kidney et", "Liver et", "Spleen et", "Gut et", "Other et",
        "Lungs", "Heart", "Kidney", "Liver", "Spleen",  "Gut", "Other"]  

#Experimental data
time_dong = [0, 5, 30, 60, 120]*60
kidney_dong_100 = [0, 1.326803, .75, .5, .25]
liver_dong_100 = [0, 64.3331574 ,80.1237696, 83.0615579 ,84.1632286]
spleen_dong_100 = [0 ,3.4820602 ,3.21158546 ,3.15035765 ,3.04021763]
heart_dong = 3.3389
lung_dong = 0.40230
# Initial NP concentration injected:
A_i=[.0019956689, .0009621096822, .0033612343, .0045212562, .0004750356, .003333333, .003333333].*1000 #mouse dm^2 
# Length of NP-cell receptor bond
l_NC_b  = (100*10^(-9)) # m

Q_i=[74.83758415, 138.673326, 763.9168803, 284.6716878, 224.48828292, 333, 1000]./(60*10^6) #L/sec MOUSE
V_i_BL=[43.65525742 ,20.04395171 ,1018.55584 ,127.7673516 ,2.226207538 ,90 ,750]./10^6
V_i_T=[33.26114851 ,16.03516137 ,56.02057122 ,75.3542703 ,7.91726 ,20 ,1000]/10^6 #L MOUSE

V_vein = 466.9/10^6 # mouse L
Q_hep = sum(Q_i[4:6]) #flow out of the liver
V_hep = sum(V_i_BL[4:6]) #total volume of blood passed thru liver
Q_vein = Q_i[3]+Q_hep+Q_i[7] #L/sec Mouse

slope1 = (22-4)/(500-50)
slope2 = 3.9
slope = slope1*slope2
a = [4, 15, 50, 79, 100]
r = a/2
K_EC0=[1.23e42, 10.6565662, 41793.3924, 2.28e12, 7.74e21, 3.33e15, 3.33e15] #162


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

K_i_on = K_i_on[5, :]
K_i_off = K_i_off[5, :]
K_deg_vein = 1/6000       

function func!(dydt, y, p, t)
  dydt[1] = (Q_i[3]*y[5]/V_i_BL[3] + Q_hep*y[6]/V_i_BL[4] + Q_i[7]*y[9]/V_i_BL[7] - Q_vein*y[1]/V_vein - y[1]*K_deg_vein) # with degradation: 
  # Artery 
  dydt[2] = (Q_vein*y[4]/V_i_BL[2]-(Q_i[3]*y[2] + Q_i[4]*y[2]+ Q_i[5]*y[2] + Q_i[6]*y[2] + Q_i[7]*y[2])/V_vein - y[2]*K_deg_vein) # 

  # Lung Vas
  dydt[3] = (Q_i[1]*y[4]/V_i_BL[2] - Q_i[1]*y[3]/V_i_BL[1]-K_i_on[1]*y[3]+K_i_off[1]*y[10]-K_NS[1]*y[3]-y[3]*K_i_deg[1])

  # Heart Vas
  dydt[4] = (Q_vein*y[1]/V_vein + Q_i[1]*y[3]/V_i_BL[1]- (Q_i[1]*y[4] + Q_vein*y[4])/V_i_BL[2] - K_i_on[2]*y[4] + K_i_off[2]*y[11]-K_NS[2]*y[4]-y[4]*K_i_deg[2])

  #Kidney Vas
  dydt[5] = (Q_i[3]*y[2]/V_vein - Q_i[3]*y[5]/V_i_BL[3] - K_i_on[3]*y[5]+K_i_off[3]*y[12]-K_NS[3]*y[5]-y[5]*K_i_deg[3])

  #Liver Vas
  dydt[6] = ((Q_i[5]*y[7]/V_i_BL[5] + Q_i[6]*y[8]/V_i_BL[6] + Q_i[4]*y[2]/V_vein - Q_hep*y[6]/V_i_BL[4] ) - K_i_on[4]*y[6] + K_i_off[4]*y[13] - K_NS[4]*y[6] - y[6]*K_i_deg[4])

  #Spleen Vas
  dydt[7] = (Q_i[5]*y[2]/V_vein - Q_i[5]*y[7]/V_i_BL[5]-K_i_on[5]*y[7]+K_i_off[5]*y[14]-K_NS[5]*y[7]-y[7]*K_i_deg[5])

  # Gut Vas
  dydt[8] = (Q_i[6]*y[2]/V_vein - Q_i[6]*y[8]/V_i_BL[6]-K_i_on[6]*y[8]+K_i_off[6]*y[15]-K_NS[6]*y[8]-y[8]*K_i_deg[6])

  # Others Vas
  dydt[9] = (Q_i[7]*y[2]/V_vein - Q_i[7]*y[9]/V_i_BL[7]-K_i_on[7]*y[9]+K_i_off[7]*y[16]-K_NS[7]*y[9]-y[9]*K_i_deg[7])
  
  # --------------- NP bound to Endothelial cell layer ---------------- #
  dydt[10:16] = (-K_i_up.*y[10:16]+K_i_on.*y[3:9]-K_i_off.*y[10:16]-y[10:16].*K_i_deg[:])

  # --------------------- NP in tissue of each organ -------------------- #
  dydt[17:23] = (K_i_up.*y[10:16] + K_NS.*y[3:9]-K_i_deg[:].*y[17:23])
  dydt[24] = (y[1]+y[2])*K_deg_vein + sum((y[3:9]+y[10:16]+y[17:23]).*K_i_deg[:])
end

y0=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0]
tspan = (0.0, 1e4)

p = [2.0e-6 1.0 1.0; 2e-6 0.1 0.02; 0.0004 0.0033 0.02; 5.0e-6 0.005 0.02; 2e-5 0.08 0.014; 3.3e-5 0.01 0.02; 3.3e-5 0.01 0.02]
#Rate of Degredation of NP in each organ
K_i_deg = p[:, 1]
#Rate of uptake of NP in each organ
K_i_up = p[:, 2]
#Rate of non specific uptake of NP in each organ
K_NS= p[:, 3]

prob = ODEProblem(func!, y0, tspan)
sol_ub = solve(prob, AutoTsit5(Rosenbrock23()))

#Plotting Mass Conservation
mass = sum(sol_ub, dims=1)
mass = vec(mass)
mass = round.(mass)
display(plot(sol_ub.t/3600, mass, ylim=(0.5, 1.5),xlabel="Time (Hours)", ylabel="Moles NP",title = "Mass Conservation",
              legend = false,size=(500,200),xtickfontsize=8,yguidefontsize=10,xguidefontsize=9, ytickfontsize=8))

#Plotiing Concentration Profiles
plot_array = Any[]
push!(plot_array, plot(0.0001.+sol_ub.t./3600, sol_ub[1:2, :]', xaxis=:log,title="Vein/Artery",legend = :right, size=(650,750), label=["Vein" "Artery"]))
for i = 17:23
  if i==17
    p1 = plot(0.0001.+sol_ub.t./3600, sol_ub[i, :]*V_vein/V_i_T[i-16], xaxis=:log)
    p2 = scatter!([float(time_dong[end])]/3600, [lung_dong], title=names_new[i],legend = false, size=(650,750))
    push!(plot_array, p2)
  elseif i==18
    p1 = plot(0.0001.+sol_ub.t./3600, sol_ub[i, :]*V_vein/V_i_T[i-16], xaxis=:log)
    p2 = scatter!([float(time_dong[end])]/3600, [heart_dong], title=names_new[i],legend = false, size=(650,750))
    push!(plot_array, p2)    
  elseif i==19
    p1 = plot(0.0001.+sol_ub.t./3600, sol_ub[i, :]*V_vein/V_i_T[i-16], xaxis=:log)
    p2 = scatter!(0.0001.+float(time_dong)/3600, kidney_dong_100, title=names_new[i],legend = false, size=(650,750))
    push!(plot_array, p2)        
  elseif i==20
    push!(plot_array, plot(0.0001.+sol_ub.t./3600, (sol_ub[i, :])*V_vein/V_i_T[i-16], xaxis=:log, title=names_new[i],legend = false, size=(650,750)))
  elseif i==21
    p1 = plot(0.0001.+sol_ub.t./3600, sol_ub[i, :]*V_vein/V_i_T[i-16], xaxis=:log)
    p2 = scatter!(0.0001.+float(time_dong)/3600, spleen_dong_100, title=names_new[i],legend = false, size=(650,750))
    push!(plot_array, p2)       
  else
    push!(plot_array, plot(0.0001.+sol_ub.t./3600, sol_ub[i, :]*V_vein/V_i_T[i-16], xaxis=:log, title=names_new[i],legend = false, size=(650,750)))
  end
end
display(plot(plot_array..., layout=(4, 2), xlabel="Time (Hours)",xtickfontsize=10, ytickfontsize=10))
