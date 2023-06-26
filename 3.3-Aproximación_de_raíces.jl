### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6bf6ed16-5315-4377-8da4-f1b4629d1c34
using Interact

# ╔═╡ c785f9de-da4f-4899-8aa7-0bb9cde22fe5
using Markdown, InteractiveUtils

# ╔═╡ 8adfc16f-60ad-4b91-a6fd-23a7232e2126
using PlutoUI

# ╔═╡ 152c306c-9aff-419c-881e-1cd6999777ce
using Pluto

# ╔═╡ 8c6bcfb6-a597-4139-afd3-d239efbb13e7
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 4a059976-8c17-11ec-20da-01784e520c36
md"""# Método de Newton

En muchas aplicaciones matemáticas, es de gran utilidad encontrar valores para los cuales una función arbitraria se anula, conocidos como las  _raíces_ de dicha función. Por ejemplo, si en un tiro parabólico podemos encontrar el valor de tiempo (positivo) en el que la función de posición vertical (altura) es igual a cero, entonces podemos determinar en qué momento la "partícula" toca el "suelo".

Más generalmente, puede ser muy útil saber para qué valores de una variable independiente $x$ una función $f$ es igual a algún valor arbitrario $b$; sin embargo, esto se puede reducir al problema anterior, pues

$$f(x) = b \quad \iff \quad f(x) - b = 0,$$

por lo que encontrar estos valores equivale a encontrar las raíces de la función $g(x) = f(x)-b$.

El _método de Newton_ (o _de Newton-Raphson_) es un método numérico iterativo capaz de encontrar raíces de una función **utilizando su derivada**, pues se sustenta en el hecho de que una función continua (la solución que buscamos para nuestro método numérico) y derivable (una suposición extra) se puede aproximar en distancias cortas como una línea recta. Este método toma como entrada un valor numérico que _sospechamos_ que podría ser una raíz -o _estar cerca_ de una raíz- de la función en cuestión y devuelve el valor aproximado de una raíz de la función.

"""

# ╔═╡ 83d0773d-4327-421c-9b54-f2d83c2a85dd
md""" ## Derivación matemática del método de Newton

Sean $f:\mathbb{R}\to\mathbb{R}$ una función derivable (y, por ende, continua), $x^\ast$ el valor exacto de una raíz de $f$ y $x_0$ una "adivinanza" de una raíz de $f$, que suponemos muy cercana al valor de $x^\ast$.

Definiendo a

$$\delta := x^\ast - x_0$$

tenemos trivialmente que

$$x^\ast = x_0 + \delta;$$

sin embargo, nuestro problema radica en que _no conocemos el valor de_ $x^\ast$ (pues, de lo contrario, no tendríamos necesidad de implementar este método numérico), por lo que no podemos calcular el valor de $\delta$. Por lo tanto, intentaremos _aproximar_ el valor de $\delta$ y, con ello, obtener una _aproximación_ de $x^\ast$.

Recordemos que queremos encontrar a $x^\ast$ tal que $f(x^\ast) = 0.$ Sustituyendo, tenemos que

$$f(x_0 + \delta) = 0.$$

Como **suponemos que $f$ es derivable** y que $x_0$ **es muy cercano a** $x^\ast$ -y que, por ende, el valor de $\delta$ es muy pequeño-, podemos expandir la ecuación anterior en una serie de Taylor alrededor de $x_0$.

"""

# ╔═╡ 00fe1fa8-855b-4287-99ca-b1839155f908
md"""

**Recordatorio** Si una función $f$ es derivable en un punto $a$ y $x$ es un valor cercano a $a$, entonces

$$f(x) = \sum_{k=0}^\infty \frac{f^k(a)(x-a)^k}{k!},$$

donde $f^k$ es la $k$-ésima derivada de $f$ y, en particular, $f^0=f$. 

Sustituyendo $x = x^\ast = x_0+\delta$ y $a=x_0$, tenemos que

$$f(x_0+\delta) = f(x_0) + f'(x_0)\delta + \dots = 0.$$

Tomando los dos primeros términos de la serie, tenemos que

$$f(x_0) + f'(x_0)\delta \approx 0,$$

de donde obtenemos la aproximación

$$\delta \approx -\frac{f(x_0)}{f'(x_0)}.$$

Por lo tanto, por definición de $\delta$, se sigue que

$$x^\ast \approx x_0 -\frac{f(x_0)}{f'(x_0)}.$$

"""

# ╔═╡ 30ff22cd-0cf3-4471-adc2-144b967e5659
md"""

Ahora que tenemos una primera aproximación de $x^\ast$ que **es más cercana a** $x^\ast$ que nuestra adivinanza inicial $x_0$, podemos aplicarle el mismo razonamiento a esta primera aproximación para encontrar una _segunda aproximación_ del valor $x^\ast$. Es decir, definimos

$$x_1 = x_0 -\frac{f(x_0)}{f'(x_0)}$$

y, aplicando los mismos argumentos, obtenemos

$$x_2 = x_1 -\frac{f(x_1)}{f'(x_1)}.$$

De aquí viene la naturaleza _iterativa_ del método de Newton.

**Ejercicio** Sean $i\geq0$ y $x_i$ la aproximación de $x^\ast$ después de $i$ iteraciones del método de Newton (la $0$-ésima iteración corresponde al valor de entrada). Encuentra una fórmula para la aproximación de $x_{i+1}$ en términos de $x_i$.

_Tu respuesta va aquí._
"""

# ╔═╡ 87922902-a3c1-4647-ab3a-eefec1721d4c
md"""## Implementación

**Ejercicio** Crea una función _newton_ que tome como argumentos a _f_, _f_$^\prime$, _x0_ y _n_, donde
* _f_ es una función,
* _f_$^\prime$ es su derivada,
* _x0_ es una aproximación inicial a una raíz de _f_, y
* _n_ es el número de iteraciones del método de Newton,

y devuelva la aproximación de una raíz de _f_ después de _n_ iteraciones.

**Sugerencia** Utiliza un ciclo iterativo _for_.

**Nota** Para obtener el nombre de función _f_$^\prime$ al escribir código de Julia, escribie _f$^\prime$_ y usa la auto completación con la tecla _<TAB>_. Esto es necesario pues _f_$^\prime$ no es un nombre válido en Julia.

"""

# ╔═╡ 3e465541-c692-4bd2-8a7f-38a883587fd2
# Tu código (comentado) va aquí :D

function newton(f, f′ , x0, n)
    x = x0

    for i in 1:n
        x -= f(x) / f′(x)
    end

    return x
end

# ╔═╡ 1c0b0804-6ff3-464e-b0a2-44b0312e2615
md"**Ejercicio** Crea una función _newtonRecursiva_ que tome los mismos argumentos que _newton_, pero implemente el método de Newton utilizando un ciclo _recursivo_.
"

# ╔═╡ 4c3d53f3-c8a3-4fd4-b549-93ae3ad5b3ed
# Tu código (comentado) va aquí :D

# Se toman los mismos parámentros de la función "newton" anterior

function newtonRecursiva(f, f′, x, n)
	
#Ahora debemos verificar si n=0
    if n == 0

# En el caso de que si sea así, se da "x"
        return x

# En caso contrario se ocupa ...
    else
        x -= f(x) / f′(x)
        return newtonRecursiva(f, f′, x, n-1)
    end
end

# ╔═╡ fa7a3432-9f19-4c80-b241-0f80c4c5dae6
md"""

**Ejercicio** Crea una función _newtonDistancia_ que tome los argumentos _f_, _f_$^\prime$, _x0_ y, en vez de _n_, tome un nuevo argumento $\varepsilon$ y realice iteraciones del método de Newton hasta que dos aproximaciones consecutivas del método tengan una distancia menor a $\varepsilon$. 

**Nota** Recordando que la distancia entre dos números reales se define como el _valor absoluto_ de la diferencia entre ellos, ¿qué sucedería si ejecutáramos una función como la anterior ingresando un valor $\varepsilon\leq0$? Por ello, haz que tu función _newtonDistancia_ imprima un mensaje de error si el valor ingresado de $\varepsilon$ es menor o igual a cero.

"""

# ╔═╡ cf003fbc-229a-4ba9-b36f-ee524c9213d9
# Tu código (comentado) va aquí :D

# Se crea la función con las mismas variables (a excepsión de (n)

function newtonDistancia(f, f′, x0, epsilon)
    if epsilon <= 0
        error("El valor de epsilon debe ser mayor que cero.")
    end
    
    x_prev = x0
    x = x_prev - f(x_prev) / f′(x_prev)
    
    while abs(x - x_prev) > epsilon
        x_prev = x
        x = x_prev - f(x_prev) / f′(x_prev)
    end
    
    return x
end

# ╔═╡ 0e62f809-0899-40e2-948c-b3d9eff1360c
md"""### Derivación numérica

La derivada de una función en un punto (suponiendo que existe) también se puede aproximar numéricamente.

Sea $f:\mathbb{R}\to\mathbb{R}$ una función derivable en un punto $x$. Entonces, la derivada de $f$ en $x$ se define como

$$f^\prime(x) = \lim_{h\to0}\frac{f(x+h)-f(x)}{h}.$$

Dado que dividir entre números muy pequeños (cercanos a cero) en una computadora puede causar errores de precisión muy grandes y $h$ es un valor totalmente arbitrario, podemos sustituir a la "variable muda" $h$ por $\frac{1}{h}$ directamente para obtener

$$f^\prime(x) = \lim_{\frac{1}{h}\to0}\frac{f(x+\frac{1}{h})-f(x)}{\frac{1}{h}}$$

o, equivalentemente,

$$f^\prime(x) = \lim_{h\to\infty} h \bigg( f\bigg(x+\frac{1}{h}\bigg)-f(x) \bigg).$$

Por lo tanto, si $h$ es un valor muy grande, tenemos que

$$f^\prime(x) \approx h \bigg(f\bigg(x+\frac{1}{h}\bigg)-f(x)\bigg). \quad (\ast)$$

"""

# ╔═╡ 30651437-7953-4527-b6aa-1fa7d2a25026
md"""**Ejercicio** Crea una función _derivadaNumérica_ que tome argumentos _f_, _x_ y _h_, donde
* _f_ es una función,
* _x_ es un punto del dominio de _f_, y
* _h_ es un valor grande,
y calcule una aproximación de $f^\prime(x)$ usando el valor `h`.

**Sugerencia** Puedes crear esta función en una sola línea (sin usar la sintáxis que utiliza la _keyword_ function).
"""

# ╔═╡ 356bdd92-50c0-4410-bd83-6c34e8951313
# Tu código (comentado) va aquí :D

derivadaNumérica(f, x, h) = (f(x + h) - f(x)) / h

# ╔═╡ 83fd7001-0fdb-4be4-b781-12bf8da68820
md"**Ejercicio** Crea una variación de la función _newton_ llamada _newtonDN_ que implemente el método de Newton calculando la **derivada numérica** de _f_ en vez de tener a la derivada de _f_ como argumento.

**Sugerencia** Usa la función _derivadaNumérica_."

# ╔═╡ 339d6090-a4a2-40d5-8812-9bf0a801ee66
# Tu código (comentado) va aquí :D

# Retomando los paràmetros de "nexton"
function newtonDN(f, x0, n, h)
    for i in 1:n
        f_prime = derivadaNumérica(f, x0, h)
        x0 -= f(x0) / f_prime
    end
    return x0
end

# ╔═╡ dcb61b75-988e-40c1-abf4-6b3092ee6be4


# ╔═╡ 5027e956-9c8c-436d-a019-b9b367ab021e
md"""**Ejercicio** Crea parámetros interactivos 
* _x0_ en el rango _-10:0.1:10_, para la "adivinanza" inicial;
* _h_ en _1000:1000:1000000_, para la aproximación de las derivadas numéricas como en la ecuación $(\ast)$;
* _n_ en _1:1000_, para el número de iteraciones del método de Newton;
y utilízalos para crear una gráfica de la función
"""

# ╔═╡ 005372a4-35cc-4d81-846e-0aa9b66a6ce5
f = cos # Luego, usa esta celda para reasignar `f`

# ╔═╡ 5e815d1d-4b27-4f22-bd33-9b8d7de2ab87
gf(x) = x^2 - 4

# ╔═╡ f8d5d895-d8c7-4a02-aa87-fccf8966efc1
md"""x0 = $(@bind x0 Slider(-10:0.1:10, 10.0, true))"""

# ╔═╡ 392502d3-db26-4870-a597-f369720eb44d
md"""h = $(@bind h Slider(1000:0.1:10000, 10000.0, true))"""

# ╔═╡ 83c53dc9-7deb-460c-ba7c-53d9b258d122
md"""n = $(@bind n Slider(1:0.1:10000, 10000.0, true))"""

# ╔═╡ 7c4e11df-fdb8-44b0-9297-2c09a6440311
begin
	x = -10:0.1:10
	y = f.(x)
	
	plot(x, y, label="f(x)", legend=:topleft)
	
	x_init = x0
	approx_raiz = newtonDN(f, x_init, n, h)
	scatter!([approx_raiz], [f(approx_raiz)], label="Aproximación de la raíz", color=:red)
end

# ╔═╡ b55c41d7-e1e1-4378-a3de-3371b07ce929
md"""

donde se muestren con una marca todos los puntos $(x_i,f(x_i))$ para cada aproximación sucesiva de la raíz, unidos por líneas rectas. ¿Qué observas para diferentes valores de $x_0$? ¿Qué sucede si le asignas funciones distintas de `cos` a la variable `f`?

**Sugerencia** Usa la función `newtonDN` y un ciclo `for`. Como cada celda de Pluto puede mostrar sólo un deslizador como salida, debes crear una celda individual por cada parámetro interactivo.

"""

# ╔═╡ 8a43221c-8c7c-485f-95fb-3f21b78cf79f
# Tu código (comentado) va aquí :D

# ╔═╡ 9a9b1dae-b120-4baa-8780-64d10f2ac0fe
md"""## Recursos complementarios

* Sección 2.3 "Newton's Method and Its Extensions" de Burden et al, _Numerical Analysis_ (2019).

"""

# ╔═╡ 2bfe89bd-cd59-4358-b0ab-e19c0ae92227
md"""## Créditos

Este _notebook_ está basado parcialmente en los _notebooks_ originales [`0.8 Raices de funciones uni-dimensionales.ipynb`](https://github.com/dpsanders/FisicaComputacional2019_3/blob/master/notebooks/08.%20Raices%20de%20funciones%20uni-dimensionales.ipynb) y [`0.9 El metodo de Newton.ipynb`](https://github.com/dpsanders/FisicaComputacional2019_3/blob/master/notebooks/09.%20El%20metodo%20de%20Newton.ipynb) del repositorio [`FisicaComputacional2019_3`](https://github.com/dpsanders/FisicaComputacional2019_3) del Dr. David Philip Sanders.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Interact = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Interact = "~0.10.5"
Pluto = "~0.19.26"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "041ade44a854af648428aa1288135cc890d9b121"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.CSSUtil]]
deps = ["Colors", "JSON", "Markdown", "Measures", "WebIO"]
git-tree-sha1 = "b9fb4b464ec10e860abe251b91d4d049934f7399"
uuid = "70588ee8-6100-5070-97c1-3cb50ed05fe8"
version = "0.1.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "62a7c76dbad02fdfdaa53608104edf760938c4ca"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.4"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.ExproniconLite]]
deps = ["Pkg", "TOML"]
git-tree-sha1 = "c2eb763acf6e13e75595e0737a07a0bec0ce2147"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.7.11"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "2613d054b0e18a3dea99ca1594e9a3960e025da4"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.7"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.Interact]]
deps = ["CSSUtil", "InteractBase", "JSON", "Knockout", "Observables", "OrderedCollections", "Reexport", "WebIO", "Widgets"]
git-tree-sha1 = "c5091992248c7134af7c90554305c600d5d9012b"
uuid = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
version = "0.10.5"

[[deps.InteractBase]]
deps = ["Base64", "CSSUtil", "Colors", "Dates", "JSExpr", "JSON", "Knockout", "Observables", "OrderedCollections", "Random", "WebIO", "Widgets"]
git-tree-sha1 = "aa5daeff326db0a9126a225b58ca04ae12f57259"
uuid = "d3863d7c-f0c8-5437-a7b4-3ae773c01009"
version = "0.10.10"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Knockout]]
deps = ["JSExpr", "JSON", "Observables", "Test", "WebIO"]
git-tree-sha1 = "91835de56d816864f1c38fb5e3fad6eb1e741271"
uuid = "bcebb21b-c2e3-54f8-a781-646b90f6d2cc"
version = "0.2.6"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "410fe4739a4b092f2ffe36fcb0dcc3ab12648ce1"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "fc8c15ca848b902015bd4a745d350f02cf791c2a"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "LoggingExtras", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "PrecompileTools", "REPL", "RegistryInstances", "RelocatableFolders", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "c4c4dac5c1332ab510e145eea59382847c51a6fb"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.26"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─8c6bcfb6-a597-4139-afd3-d239efbb13e7
# ╟─4a059976-8c17-11ec-20da-01784e520c36
# ╟─83d0773d-4327-421c-9b54-f2d83c2a85dd
# ╟─00fe1fa8-855b-4287-99ca-b1839155f908
# ╟─30ff22cd-0cf3-4471-adc2-144b967e5659
# ╠═87922902-a3c1-4647-ab3a-eefec1721d4c
# ╠═3e465541-c692-4bd2-8a7f-38a883587fd2
# ╠═1c0b0804-6ff3-464e-b0a2-44b0312e2615
# ╠═4c3d53f3-c8a3-4fd4-b549-93ae3ad5b3ed
# ╠═fa7a3432-9f19-4c80-b241-0f80c4c5dae6
# ╠═cf003fbc-229a-4ba9-b36f-ee524c9213d9
# ╠═0e62f809-0899-40e2-948c-b3d9eff1360c
# ╠═30651437-7953-4527-b6aa-1fa7d2a25026
# ╠═356bdd92-50c0-4410-bd83-6c34e8951313
# ╠═83fd7001-0fdb-4be4-b781-12bf8da68820
# ╠═339d6090-a4a2-40d5-8812-9bf0a801ee66
# ╠═dcb61b75-988e-40c1-abf4-6b3092ee6be4
# ╠═5027e956-9c8c-436d-a019-b9b367ab021e
# ╠═005372a4-35cc-4d81-846e-0aa9b66a6ce5
# ╠═6bf6ed16-5315-4377-8da4-f1b4629d1c34
# ╠═c785f9de-da4f-4899-8aa7-0bb9cde22fe5
# ╠═8adfc16f-60ad-4b91-a6fd-23a7232e2126
# ╠═152c306c-9aff-419c-881e-1cd6999777ce
# ╠═5e815d1d-4b27-4f22-bd33-9b8d7de2ab87
# ╠═f8d5d895-d8c7-4a02-aa87-fccf8966efc1
# ╠═392502d3-db26-4870-a597-f369720eb44d
# ╠═83c53dc9-7deb-460c-ba7c-53d9b258d122
# ╠═7c4e11df-fdb8-44b0-9297-2c09a6440311
# ╟─b55c41d7-e1e1-4378-a3de-3371b07ce929
# ╠═8a43221c-8c7c-485f-95fb-3f21b78cf79f
# ╟─9a9b1dae-b120-4baa-8780-64d10f2ac0fe
# ╟─2bfe89bd-cd59-4358-b0ab-e19c0ae92227
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
