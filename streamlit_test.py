import streamlit as st

# theory
# https://www.youtube.com/watch?v=fWgrHVAroto&list=PLWVKUEZ25V94kdT2Lh97KqB9MoLV9ZzmU&index=8

# code
# https://www.youtube.com/watch?v=BJPeXGyj6V8&list=PLWVKUEZ25V94kdT2Lh97KqB9MoLV9ZzmU&index=9


import matplotlib
# matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np


# st.set_page_config(layout="centered")
css='''
<style>
    section.main > div {max-width:90rem}
</style>
'''
st.markdown(css, unsafe_allow_html=True)

st.title("Three gene model oscillations.")
st.write("Model explanation: Gene1 activates Gene2, Gene2 activates Gene3 & Gene3 represses Gene1")
st.image("4.3 gene-model-network-oscillating gene network.jpg",
         caption="G1 activates G2, G2 activates G3, G3 represses G1",
         width=700)


st.sidebar.write("**Play with the parameters and take a look at the behaviour of the model.**")

# Initial Conditions Block
st.sidebar.write("Initial Conditions:")
col1, col2, col3 = st.columns([1,1,1])
G1 = st.sidebar.number_input("Gen1 Population", min_value=0,  value=0,  step=1)
G2 = st.sidebar.number_input("Gen2 Population", min_value=0,  value=0,  step=1)
G3 = st.sidebar.number_input("Gen3 Population", min_value=0,  value=0,  step=1)
y0 = [G1,G2,G3]
# Time
st.sidebar.write("Time:")
col1, col2, col3 = st.columns(3)
start = st.sidebar.number_input("Starting point", min_value=0,  value=0,  step=1)
end = st.sidebar.number_input("Ending point", min_value=0,  value=200, step=1)
points = st.sidebar.number_input("Number of points", min_value=0,  value=100,  step=1)

t = np.linspace(start, end, points)

st.sidebar.write("Parameters:")
k_1 = st.sidebar.number_input("Gen1 Production rate (k1)", min_value=0.0,  value=round(0.5,1), format="%1f", step=0.1)
k_2 = st.sidebar.number_input("Gen2 Production rate (k2)", min_value=0.0,  value=0.5, format="%1f", step=0.1)
k_3 = st.sidebar.number_input("Gen3 Production rate (k3)", min_value=0.0,  value=0.5, format="%1f", step=0.1)
gamma_1 = st.sidebar.number_input("Gen1 degradation rate (γ1)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
gamma_2 = st.sidebar.number_input("Gen2 degradation rate(γ2)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
gamma_3 = st.sidebar.number_input("Gen3 degradation rate (γ3)", min_value=0.0,  value=0.1, format="%1f", step=0.1)
n = st.sidebar.number_input("constant", min_value=0.0,  value=9.0, format="%1f", step=0.1)
c = st.sidebar.number_input("power", min_value=0.0,  value=1.0, format="%1f", step=0.1)


params = [k_1, gamma_1, k_2, gamma_2, k_3, gamma_3, n, c]

def sim(variables, t, params):

    G1 = variables[0]
    G2 = variables[1]
    G3 = variables[2]

    k_1 = params[0]
    gamma_1 = params[1]
    k_2 = params[2]
    gamma_2 = params[3]
    k_3 = params[4]
    gamma_3 = params[5]
    n = params[6]
    c = params[7]

    dG1dt = ((c**n) / (c**n + G3**n)) * k_1 - (gamma_1 *G1)
    dG2dt = ((G1**n) / (c**n + G1**n)) * k_2 - (gamma_2 *G2)
    dG3dt = ((G2**n) / (c**n + G2**n)) * k_3 - (gamma_3 *G3)

    return([dG1dt,dG2dt,dG3dt])

y = odeint(sim, y0, t, args=(params,))

# plotting
f, (ax1,ax2,ax3) = plt.subplots(3, sharex=True, sharey=False)

line1, = ax1.plot(t, y[:,0], color="b", label="G1")
line2, = ax2.plot(t, y[:,1], color="r", label="G2")
line3, = ax3.plot(t, y[:,2], color="g", label="G3")

ax1.set_ylabel("Number")
ax1.set_xlabel("Time")
ax1.legend(handles = [line1,line2,line3])

st.pyplot(f)
