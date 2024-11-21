# theory and code 
# https://www.youtube.com/watch?v=pzkutGeVYlM&list=PLWVKUEZ25V94kdT2Lh97KqB9MoLV9ZzmU&index=13

import matplotlib
# matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import random
import streamlit as st

# st.set_page_config(layout="centered")
css='''
<style>
    section.main > div {max-width:90rem}
</style>
'''
st.markdown(css, unsafe_allow_html=True)

st.title("Three genes model dynamics with gillespie algorithm.")
st.write("Model explanation: Gene1 activates Gene2, Gene2 activates Gene3 & Gene3 represses Gene1")
st.image("images/6. Three gene model network (stochastic).png",
         caption="Model depends on the rates of events",
         width=700)

st.sidebar.write("**Play with the parameters and take a look at the behaviour of the model.**")
# Initial Conditions
st.sidebar.write("Initial Conditions: g1=0, G2=0, G3=0, time=0")

G1 = [0] # gene1 population
G2 = [0] # gene2 population
G3 = [0] # gene3 population
t = [0] # time

#parameters
st.sidebar.write("Parameters:")
t_end = st.sidebar.number_input("end time of the simulation", min_value = 1, value = 1000, step=1)
k_1 = st.sidebar.number_input("Gen1 Production rate (k1)", min_value=0.0,  value=round(0.5,1), format="%1f", step=0.1)
k_2 = st.sidebar.number_input("Gen2 Production rate (k2)", min_value=0.0,  value=0.5, format="%1f", step=0.1)
k_3 = st.sidebar.number_input("Gen3 Production rate (k3)", min_value=0.0,  value=0.5, format="%1f", step=0.1)
gamma_1 = st.sidebar.number_input("Gen1 degradation rate (γ1)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
gamma_2 = st.sidebar.number_input("Gen2 degradation rate(γ2)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
gamma_3 = st.sidebar.number_input("Gen3 degradation rate (γ3)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
n = st.sidebar.number_input("constant", min_value=0.0,  value=9.0, format="%1f", step=0.1)
c = st.sidebar.number_input("power", min_value=0.0,  value=1.0, format="%1f", step=0.1)

while t[-1] < t_end:

    current_G1 = G1[-1]
    current_G2 = G2[-1]    
    current_G3 = G3[-1]

    rates = [(c**n / (c**n + current_G3**n)) * k_1, gamma_1 * current_G1, \
    (current_G1**n / (c**n + current_G1**n)) * k_2, gamma_2 * current_G2, \
    (current_G2**n / (c**n + current_G2**n)) * k_3, gamma_3 * current_G3, \
    ]

    rate_sum = sum(rates)

    tau = np.random.exponential(scale=1/rate_sum) # mean = 1/λ λ=sum(rates)
    t.append(t[-1] + tau) # add random T to the last timepoint

    rand = random.uniform(0,1)

    # probability of events

     # G1 production event
    if rand * rate_sum <= rates[0]:
            G1.append(G1[-1] + 1)
            G2.append(G2[-1])
            G3.append(G3[-1])

    # G1 decay event
    elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
            G1.append(G1[-1] - 1)
            G2.append(G2[-1])
            G3.append(G3[-1])

    # G2 production event
    elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
            G1.append(G1[-1])
            G2.append(G2[-1] + 1)
            G3.append(G3[-1])

    # G2 decay event
    elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
            G1.append(G1[-1])
            G2.append(G2[-1] - 1)
            G3.append(G3[-1])

    # G3 production event
    elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
            G1.append(G1[-1])
            G2.append(G2[-1])
            G3.append(G3[-1] + 1)

    # G3 decay event
    elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
            G1.append(G1[-1])
            G2.append(G2[-1])
            G3.append(G3[-1] - 1)



f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
line1, = ax1.plot(t , G1, color="b",label="G1")
line2, = ax2.plot(t , G2, color="r",label="G2")
line3, = ax3.plot(t , G3, color="g",label="G3")
ax1.set_ylabel('Number')
ax1.set_xlabel('Time')
ax2.legend(handles=[line1,line2,line3],bbox_to_anchor=(1., 0.5), loc='center left')
f.tight_layout()
st.pyplot(f)
# ax1.legend(bbox_to_anchor=(1., 0.5), loc='center left')
# plt.show()

