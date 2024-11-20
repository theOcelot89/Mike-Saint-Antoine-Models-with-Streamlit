import streamlit as st

# theory
# https://www.youtube.com/watch?v=NnLxfBWfJEE&list=PLWVKUEZ25V94kdT2Lh97KqB9MoLV9ZzmU&index=11&pp=iAQB

# code
# https://www.youtube.com/watch?v=Hgjh3YxC01o&list=PLWVKUEZ25V94kdT2Lh97KqB9MoLV9ZzmU&index=12

import matplotlib
# matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import random

css='''
<style>
    section.main > div {max-width:90rem}
</style>
'''
st.markdown(css, unsafe_allow_html=True)

st.title("Gillespie algorithm on mRNA lifecycle.")
st.write("mRNA production depends on k and degradation on γX")
st.image("images/5. gillespie algorithm.png",
         caption="Above: differential equations for deterministic model (not used), below: events and rates",
         width=700)

st.sidebar.write("**Play with the parameters and take a look at the behaviour of the model.**")

# Initial Conditions
st.sidebar.write("Initial Conditions: mRNA=0, time=0")
x = [0] # mRNA population
t = [0] # time

st.sidebar.write("Parameters:")
t_end = st.sidebar.number_input("end time of simulation", min_value = 1, value = 1000, step=1)
k = st.sidebar.number_input("mRNA production rate", min_value=0.0,  value=round(2.0,1), format="%1f", step=0.1)
gamma = st.sidebar.number_input("Decay rate", min_value=0.0,  value=round(0.1,1), format="%1f", step=0.1)



while t[-1] < t_end: # while not end of time

    current_x = x[-1]
    rates = [k, gamma * current_x]
    rates_sum = sum(rates)

    tau = np.random.exponential(scale=1/rates_sum) # mean = 1/λ λ=sum(rates)
    t.append(t[-1] + tau) # add random T to the last timepoint

    rand = random.uniform(0,1)

    # probability of events

    # production event
    if rand * rates_sum > 0 and rand * rates_sum < rates[0]: # bigger than zero and less than k
        x.append(x[-1] + 1)



    # decay event    
    elif rand * rates_sum > rates[0] and rand * rates_sum < rates[0] + rates[1]: # bigger than k and less tha γX
        x.append(x[-1] - 1)

print(sum(x)/len(x)) # average mRNA quantity, is always close to steady state  (x when nothing changes) of model  (k/γ)
st.sidebar.write(f"Mean quantiny of mRNA is **{round(sum(x)/len(x),2)}** \
                    which is close to the steady state of the deterministic version\
                    of the model (k/γ). In the deterministic model if x is \
                    **{round(sum(x)/len(x))}** then no change will occur.")

f = plt.figure()
plt.plot(t,x)
plt.xlabel("time")
plt.ylabel("mRNA quantity")
st.pyplot(f)