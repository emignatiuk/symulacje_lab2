import random
import numpy as np
import json

with open("model.json") as f:
    model = json.load(f)
    
def stworz_funkcje_szybkosci(reakcja):
    zmienne = reakcja["params"] + reakcja["species"]
    tekst_lambdy = f"lambda {', '.join(zmienne)}: {reakcja['equation']}"
    funkcja = eval(tekst_lambdy)
    return funkcja, zmienne

def symulacja_gillespie(model, max_kroki, max_czas=None):
    gatunki = list(model['species_init'].keys())
    parametry = model['parameters']
    indeksy = {nazwa: i for i, nazwa in enumerate(gatunki)}
    liczba_gatunkow = len(gatunki)

    trajektoria_gatunkow = np.zeros((max_kroki + 1, liczba_gatunkow))
    for i, nazwa in enumerate(gatunki):
        trajektoria_gatunkow[0, i] = model['species_init'][nazwa]
    trajektoria_czasu = np.zeros(max_kroki + 1)

    reakcje = []
    for rek in model['reactions'].values():
        funkcja, zmienne = stworz_funkcje_szybkosci(rek)
        reakcje.append({
            'funkcja': funkcja,
            'zmienne': zmienne,
            'efekty': rek['effects']
        })

    wykonane_kroki = 0

    for krok in range(max_kroki):
        aktualny_czas = trajektoria_czasu[krok]
        if max_czas is not None and aktualny_czas >= max_czas:
            break

        aktualny_stan = trajektoria_gatunkow[krok]
        szybkosci = []

        for rek in reakcje:
            wartosci = []
            for zm in rek['zmienne']:
                if zm in parametry:
                    wartosci.append(parametry[zm])
                else:
                    wartosci.append(aktualny_stan[indeksy[zm]])
            try:
                r = rek['funkcja'](*wartosci)
            except:
                r = 0
            szybkosci.append(max(r, 0))

        suma_szybkosci = sum(szybkosci)
        if suma_szybkosci == 0:
            break

        dt = random.expovariate(suma_szybkosci)
        trajektoria_czasu[krok + 1] = aktualny_czas + dt
        trajektoria_gatunkow[krok + 1] = aktualny_stan.copy()

        prog = random.uniform(0, suma_szybkosci)
        akumulacja = 0
        wybrana = -1
        for i, r in enumerate(szybkosci):
            akumulacja += r
            if prog <= akumulacja:
                wybrana = i
                break

        for sp, zmiana in reakcje[wybrana]['efekty'].items():
            idx = indeksy[sp]
            trajektoria_gatunkow[krok + 1][idx] += zmiana

        wykonane_kroki = krok + 1

    print(f"Wykonano kroków: {wykonane_kroki}")
    return trajektoria_czasu[:wykonane_kroki + 1], trajektoria_gatunkow[:wykonane_kroki + 1], gatunki
