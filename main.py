import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# dane ze Stellarium (Układ Równikowy, pomiary poziome przez 100 godzin
# czas od rozpoczęcia pomiaru
czas_h = np.array([
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
    64, 68, 72, 76, 80, 84, 88, 92, 96, 100
])

# odległość kątowa w sekundach łuku mierzona tylko w poziomie
# Znak '-' to prawa strona od Jowisza, '+' to lewa strona. To czy lewa czy prawa storna ma miuns nie ma dużego znacznie.
# po porstu musimy się trzymać jednego założenia
kat_arcsec = np.array([
    2.10, 44.70, 84.20, 116.70, 139.20, 149.7, 147.30, 132.10, 105.30, 69.20, 27.10, -17.6, -60.8, -98.4,
    -127.4, -145.4, -150.9, -143.6, -124.1, -94.3, -56.6, -14.2, 29.5, 70.6, 105.7, 131.9
])

# zamiana odległóśći pomiedzy zmienia a Jowiszem. wartość odczytana ze Stellarium
D_metrow = 913355000 * 1000


# zamiana (układ kątowy na układ SI)
t_dane = czas_h * 3600                   #godziny na sekundy
kat_radiany = np.radians(kat_arcsec / 3600)
x_dane = D_metrow * np.tan(kat_radiany)  # Rzeczywista odległość Europy w metrach


# parametr 'C' (offset) - na wypadek, gdyby środek Jowisza lekko drgał
def sinusoida_model(t, R, T, phi, C):
    return R * np.sin(2 * np.pi * t / T + phi) + C

# szacowane parametry startowe dla Europy
# R_start = 670 000 km, T_start = 3.55 dnia (306720 sekund), phi_start = 0, C_start = 0
szacunki_startowe = [670000000, 306720, 0, 0]

# dopasowanie krzywej metodą najmniejszych kwadratów
parametry, _ = curve_fit(sinusoida_model, t_dane, x_dane, p0=szacunki_startowe)
R_wyznaczone, T_wyznaczone, phi_wyznaczone, C_wyznaczone = parametry


R_wyznaczone = abs(R_wyznaczone)
T_wyznaczone = abs(T_wyznaczone)

G = 6.6743e-11
M_Jowisza = (4 * np.pi**2 * R_wyznaczone**3) / (G * T_wyznaczone**2) # III prao Keplra


plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(11, 6.5))


t_gładkie = np.linspace(t_dane.min(), t_dane.max(), 1000)
x_gładkie = sinusoida_model(t_gładkie, R_wyznaczone, T_wyznaczone, phi_wyznaczone, C_wyznaczone)


ax.scatter(czas_h, x_dane / 1000, color='#00ffcc', label='Pomiary ze Stellarium (poziome)', s=50, zorder=3)
ax.plot(t_gładkie / 3600, x_gładkie / 1000, color='#ff007f', label='Teoretyczny tor ruchu (Dopasowana sinusoida)', linewidth=2)


ax.axhline(C_wyznaczone / 1000, color='white', linestyle='--', alpha=0.3, label='Środek układu (Jowisz)')
ax.set_title('Projekt: Wyznaczanie masy Jowisza z III Prawa Keplera', fontsize=14, color='yellow', pad=20)
ax.set_xlabel('Czas od rozpoczęcia serii pomiarowej (h)', fontsize=11)
ax.set_ylabel('Pozycja Europy względem Jowisza (km)', fontsize=11)
ax.grid(True, linestyle=':', alpha=0.3)
ax.legend(loc='upper right', framealpha=0.8)


tekst_raportu = (
    f"📊 WYNIKI ANALIZY:\n"
    f"─────────────────────────────────────\n"
    f"Promień orbity (R):  {R_wyznaczone/1000:.1f} km\n"
    f"Okres obiegu (T):    {T_wyznaczone/3600:.2f} h ({T_wyznaczone/(3600*24):.3f} dni)\n"
    f"Odchyłka środka (C): {C_wyznaczone/1000:.2f} km\n"
    f"─────────────────────────────────────\n"
    f"🔥 WYLICZONA MASA JOWISZA:\n"
    f"M = {M_Jowisza:.4e} kg\n\n"
    f"Wartość tablicowa: 1.8983e+27 kg"
)


plt.gcf().text(0.14, 0.15, tekst_raportu, fontsize=10, fontfamily='monospace',
               bbox=dict(facecolor='#111111', alpha=0.9, edgecolor='yellow', boxstyle='round,pad=1'))

plt.tight_layout()
plt.show()