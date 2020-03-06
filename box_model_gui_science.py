# -------------------------------------------------------------------
# Box model of nutrient cycling the Lower St Lawrence Estuary (LSLE)
#
# Designed at McGill University, copyright 2020
# Mathilde Jutras
#
# The LSLE is a three-layer stratified system.
# This GUI program finds the concentration of nutrients in each layer
# under perturbations in the input of nutrient.
#
# For details on the model, refer to Jutras, M. et al., (under review),
# Nutrient cycling in the Lower St. Lawrence Estuary: response to 
# environmental perturbations, Estarine, Coastal and Shelf Science.
# -------------------------------------------------------------------

import Tkinter as tk
import tkMessageBox
import numpy as np

# ------------------------------------------
# FUNCTIONS
# -----------------------------------------

# Get the particulate flux necessary to close the nutrient budget
def get_Fp( cs1_d, cs1_p, cI, cD, c1, S, fs1, fI, fD, fs2 ) :

	# Output flux at the surface towards the Gulf
	Fp = fs1*cs1_d + fs1*cs1_p + fI*cI + fD*cD - fs2*c1 - S
	# Concentration of the output
	c1p = Fp / fs2

	return Fp, c1p

# Get the settling rates
def settling( cs1_d, cs1_p, c1, c1p, cD, c3, c2, S, fs1, fI, fD, fs2 ) :

	fs2 = fs1+fI+fD 	# output at surface
	fDI = fD 			# upwelling deep to intermediate
	fIS = fI+fD 		# upwelling intermedate to surface

	# Solve the system of equations
	a = np.array([[1, 0], [0, -1]])
	b = np.array([fs1*(cs1_d+cs1_p) - fs2*(c1+c1p) + fIS*c2, fD*cD - fDI*c3 - S])
	P1, P2 = np.linalg.solve(a,b)

	# Get this in terms of the particulate content of the first layer (P = a*c1,p)
	a1 = P1 / c1p
	a2 = P2 / c1p

	return P1, P2, a1, a2

# Calculate the mineralisation and remineralisation rates
def mineral( fs1, cs1_p, cs1_d, fs2, fI, fD, cI, cD, c1p, c1, c2, c3, P1, P2, S ) :

	fs2 = fs1+fI+fD 	# output at surface
	fDI = fD			# upwelling deep to intermediate
	fIS = fI+fD 		# upwelling intermedate to surface

	# Nutrient uptake at the surface. Ga should equal -Gb
	Ga = fs1*cs1_p - fs2*c1p - P1
	Gb = fs1*cs1_d - fs2*c1 + fIS*c2
	alpha = abs(Ga/c1)  # G = alpha * C1d
	# Remineralisation rates in intermediate (R1) and deep (R2) layers
	# R1a should equal -R1b and R2a should equal -R2b
	R1a = P1-P2
	R1b = fDI*c3 - fIS*c2 + fI*cI
	R2a = P2-S
	R2b = fD*cD - fDI*c3

	return Ga, Gb, alpha, R1a, R1b, R2a, R2b

# Solve the model to find nutrient concentrations in each layer
def solve_model( cs1p, cs1d, cI, cD, fs1, fI, fD, a1, a2, bS, alpha ) :

	fs2=fs1+fI+fD

	a = np.array([[fs2, fs2+a1, -fI-fD, 0], [0, -a2+bS, 0, fD], [fs2, fs2+bS, 0, 0], [-alpha, fs2+a1, 0, 0]])
	b = np.array([fs1*(cs1p+cs1d), fD*cD, fs1*(cs1p+cs1d)+fI*cI+fD*cD, fs1*cs1p])

	# c1d = surface, dissolved
	# c1p = surface, particulate
	# c2  = intermediate, dissolved
	# c3  = deep, dissolves

	c1d, c1p, c2, c3 = np.linalg.solve(a,b)

	return c1p, c1d, c2, c3

# Main script
def run_model() :

	# Check that have numbers
	try:
		float(Cs1_d.get()) ; float(Cs1_p.get()) ; float(CI.get()) ; float(CD.get())
	except ValueError:
		tkMessageBox.showinfo('Window Title', u'Les valeurs entrées devraient être des nombres')

	# --------------------------------------
	# INITIAL CONDITIONS
	# --------------------------------------

	# Volume each layer
	V1 = 40.*1000. * 50. * 200. * 1000.
	V2 = 15.*1000. * 100. * 200. * 1000.
	V3 = 15.*1000. * 150. * 200. * 1000.

	# Input volume flux
	fs1 = 1.39*10**4 	# surface input
	fD = 3.*10**4 		# deep layer input
	fI = 4.*10**4 		# intermediary layer input
	fs2 = fs1+fI+fD 	# output at surface

	# Nitrate
	cN1 = 10.*10**-3  	# surface layer
	cN2 = 15.*10**-3	# intermediate layer
	cN3 = 23.*10**-3	# deep layer
	cNs1_d = 23.4*10**-3 # surface dissolved input
	cNs1_p = 22.5*10**-3 # surface particulate input
	cNI = 7.23*10**-3 	# intermediate input
	cND = 21.5*10**-3 	# deep input

	# Phosphate
	cP1 = 0.82*10**-3
	cP2 = 1.33*10**-3
	cP3 = 1.72*10**-3
	cPs1_d = 4.52*10**-4
	cPs1_p = 5.81*10**-4
	cPD = 1.50*10**-3
	cPI = 0.93*10**-3

	# Silice
	cSi1 = 14.*10**-4
	cSi2 = 18.*10**-3
	cSi3 = 35.*10**-3
	cSis1_d = 42.6*10**-3
	cSis1_p = 11./21. * cSis1_d
	cSiI = 5.69*10**-3
	cSiD = 14.5*10**-3

	# Permanent burial rates
	SN = 1.4
	SP = 3.8
	SSi = 119.

	# --- Particulate flux out --- #
	# Obtain the particulate surface output necessary to ensure nutrient budget closing

	[FpN, cN1p]   = get_Fp( cNs1_d, cNs1_p, cNI, cND, cN1, SN, fs1, fI, fD, fs2 )
	[FpP, cP1p]   = get_Fp( cPs1_d, cPs1_p, cPI, cPD, cP1, SP, fs1, fI, fD, fs2 )
	[FpSi, cSi1p] = get_Fp( cSis1_d, cSis1_p, cSiI, cSiD, cSi1, SSi, fs1, fI, fD, fs2 )

	# Settling rates
	P1N, P2N, a1N, a2N = settling( cNs1_d, cNs1_p, cN1, cN1p, cND, cN3, cN2, SN, fs1, fI, fD, fs2 )
	P1P, P2P, a1P, a2P = settling( cPs1_d, cPs1_p, cP1, cP1p, cPD, cP3, cP2, SP, fs1, fI, fD, fs2 )
	P1Si, P2Si, a1Si, a2Si = settling( cSis1_d, cSis1_p, cSi1, cSi1p, cSiD, cSi3, cSi2, SSi, fs1, fI, fD, fs2 )

	# Mineralisation/remineralisation rates
	[GNa, GNb, alphaN, RINa, RINb, RDNa, RDNb] = mineral( fs1, cNs1_p, cNs1_d, fs2, fI, fD, cNI, cND, cN1p, cN1, cN2, cN3, P1N, P2N, SN )
	[GPa, GPb, alphaP, RIPa, RIPb, RDPa, RDPb] = mineral( fs1, cPs1_p, cPs1_d, fs2, fI, fD, cPI, cPD, cP1p, cP1, cP2, cP3, P1P, P2P, SP )
	[GSia, GSib, alphaSi, RISia, RISib, RDSia, RDSib] = mineral( fs1, cSis1_p, cSis1_d, fs2, fI, fD, cSiI, cSiD, cSi1p, cSi1, cSi2, cSi3, P1Si, P2Si, SSi )

	# Express burial rates in term of particulate concentration in the first layer
	bSN = SN / cN1p
	bSP = SP / cP1p
	bSSi = SSi / cSi1p

	# --------------------------------------
	# CALCULATE NEW STEADY-STATE
	# --------------------------------------

	# Get the values from the GUI
	cNs1d = float(Cs1_d.get()) ; cNs1p = float(Cs1_p.get()) ; cNI = float(CI.get()) ; cND = float(CD.get())

	[FpN, cN1p]   = get_Fp( cNs1d, cNs1p, cNI, cND, cN1, SN, fs1, fI, fD, fs2 )

	# Solve the model
	c1pN, c1dN, c2N, c3N = solve_model(cNs1p, cNs1d, cNI, cND, fs1, fI, fD, a1N, a2N, bSN, alphaN)
	P1N = a1N*c1pN ; P2N = a2N*c1pN
	[GaN, GbN, alphaN, RIaN, RIbN, RDaN, RDbN] = mineral( fs1, cNs1p, cNs1d, fs2, fI, fD, cNI, cND, cN1p, cN1, cN2, cN3, P1N, P2N, SN )

	#c1pP, c1dP, c2P, c3P = solve_model(cPs1p, cPs1d, cPI, cPD, fs1, fI, fD, a1P, a2P, bSP, alphaP)
	#P1P = a1P*c1pP ; P2P = a2P*c1pP
	#[GaP, GbP, alphaP, RIaP, RIbP, RDaP, RDbP] = mineral( fs1, cPs1p, cPs1d, fs2, fI, fD, cPI, cPD, c1pP, c1dP, c2P, c3P, P1P, P2P, SP )

	#O2cons3 = min( 150.*RDaP, 150./16.*RDaN )
	O2cons3 = 150./16.*RDaN
	O2cons3 = O2cons3/V3*3600.*24.*365.25

	# Inclue effect of a change in temperature
	try:
		float(dT.get())
		if float(dT.get()) > 0. :
			O2cons3 = O2cons3 * 2.55**(float(dT.get())/10.)
	except ValueError:
		O2cons3 = O2cons3

# --------------------------------------------------
# SCRIPT FOR GUI
# --------------------------------------------------

root = tk.Tk()
frame = tk.Frame(root)
frame.pack()

l = 0

instr = tk.Label(frame, text='Eutrophication of the St. Lawrence Estuary')
instr.grid(row=l, columnspan=3) ; l+=2
instr.config(font=('Arial',14, 'bold','underline'))

instr = tk.Label(frame, text=' ') # empty line
instr.grid(row=l, columnspan=3) ; l+=2

instr = tk.Label(frame, text='Increases in the amount of nitrogen-based nutrients reaching the St. Lawrence Estuary foster eutrophication. \n This program evaluates the approximate impact of changing nitrate concentration in the different inputs to \n the system.', justify='left')
instr.grid(row=l, columnspan=3, sticky='w') ; instr.config(font=('Arial',12)) ; l+=2

instr = tk.Label(frame, text=' ') # empty line
instr.grid(row=l, columnspan=3) ; l+=2

t1 = tk.Label(frame, text='Current value') ; t1.grid(row=l, column=1) ; t1.config(font=('Arial',12,'bold'))
t2 = tk.Label(frame, text='New value') ; t2.grid(row=l, column=2) ; t2.config(font=('Arial',12,'bold'))
l+=1
t1_1 = tk.Label(frame, text='(micromol/L)') ; t1_1.grid(row=l, column=1) ; t1_1.config(font=('Arial',12,'bold'))
t2_1 = tk.Label(frame, text='(micromol/L)') ; t2_1.grid(row=l, column=2) ; t2_1.config(font=('Arial',12,'bold'))
l+=1

# Get data and results

l1 = tk.Label(frame, text='St Lawrence River input - dissolved concentration')
l1.grid(row=l+1,column=0, sticky='w')
l1.config(font=('Arial', 12))
l2 = tk.Label(frame, text='St Lawrence River input - particulate concentration')
l2.grid(row=l+2,column=0, sticky='w')
l2.config(font=('Arial', 12))
l3 = tk.Label(frame, text='Intermediary input concentration (from Atlantic)')
l3.grid(row=l+3,column=0, sticky='w')
l3.config(font=('Arial', 12))
l4 = tk.Label(frame, text='Deep input concentration (from Atlantic)')
l4.grid(row=l+4,column=0, sticky='w')
l4.config(font=('Arial', 12))

l1_1 = tk.Label(frame, text='23.4') ; l1_1.grid(row=l+1,column=1) ; l1_1.config(font=('Arial',12))
l2_1 = tk.Label(frame, text='22.5') ; l2_1.grid(row=l+2,column=1) ; l2_1.config(font=('Arial',12))
l3_1 = tk.Label(frame, text='7.4') ; l3_1.grid(row=l+3,column=1) ; l3_1.config(font=('Arial',12))
l4_1 = tk.Label(frame, text='22.1') ; l4_1.grid(row=l+4,column=1) ; l4_1.config(font=('Arial',12))

Cs1_d = tk.Entry(frame)
Cs1_p = tk.Entry(frame)
CI = tk.Entry(frame)
CD = tk.Entry(frame)

Cs1_d.grid(row=l+1, column=2)
Cs1_p.grid(row=l+2, column=2)
CI.grid(row=l+3, column=2)
CD.grid(row=l+4, column=2)

l+=5

# Temperature increase
instr = tk.Label(frame, text=' ') # empty line
instr.grid(row=l, columnspan=3) ; l+=1

T1 = tk.Label(frame, text='If you wish to include the effect of a temperature increase in the deep water, do se here.')
T1.grid(row=l,column=0, columnspan=3,sticky='w') ; T1.config(font=('Arial', 12)) ; l+=1

T2 = tk.Label(frame, text='Temperature increase:')
T2.grid(row=l,column=0, sticky='e') ; T2.config(font=('Arial', 12))
dT = tk.Entry(frame) ; dT.grid(row=l, column=1)
l+=1

# Buttons
button = tk.Button(frame, text='RUN', command=run_model).grid(row=l,column=1,sticky=tk.W,pady=4)
button = tk.Button(frame, text='QUIT', command=frame.quit).grid(row=l,column=2,sticky=tk.W,pady=4)
l+=1

instr = tk.Label(frame, text=' ') # empty line
instr.grid(row=l, columnspan=3) ; l+=1

# Consumption rate result
r1 = tk.Label(frame, text='Initial oxygen consumption rate at depth:')
r1.grid(row=l, column=0, columnspan=2, sticky='w') ; r1.config(font=('Arial',12))
r2 = tk.Label(frame, text='New oxygen consumption rate:')
r2.grid(row=l+1, column=0, sticky='w') ; r2.config(font=('Arial',12,'bold'))
r1_2 = tk.Label(frame, text='40.0') ; r1_2.grid(row=l, column=1) ; r1_2.config(font=('Arial',12))
r1_3 = tk.Label(frame, text='micromol/L/yr') ; r1_3.grid(row=l, column=2) ; r1_3.config(font=('Arial',12))
r2_3 = tk.Label(frame, text='micromol/L/yr') ; r2_3.grid(row=l+1, column=2) ; r2_3.config(font=('Arial',12))

res = tk.Label(frame, bg='black', fg='white', width=8)
res.grid(row=l+1, column=1) ; res.config(font=('Arial',12,'bold'))

l+=2

t4 = tk.Label(frame, text='Values above are representative of the deep waters in front of Rimouski.') ; t4.grid(row=l,column=0,columnspan=3, sticky='w') ; t4.config(font=('Arial', 12)) ; l+=1

cr = tk.Label(frame, text=u"\u00A9" + 'Mathilde Jutras - McGill University')
cr.grid(row=l, column=2) ; cr.config(font=('Arial',10))

root.mainloop()
