### Modelo garch de volatiidad dde COP/USD ###

rm(list = ls())

pacman::p_load(
  rugarch, fGarch, aTSA,
  quantmod, forecast,
  tidyverse, tseries,     
  urca, FinTS, MASS, TSA, car 
)

# ─── 1. Importación de datos ──────────────────────────────────────────────── #

#---Para cambios de corte temporal, tanto de apertura y de cierre, cambios para las variables
#--- From y To, añadiendo la fecha en formato Americano, "Año-Mes-Día

USDCOP <- getSymbols("COP=X",
                     from        = "2025-07-01",
                     periodicity = "daily",
                     source      = "yahoo",
                     auto.assign = FALSE)[, 4]   # Precio de cierre

USDCOP <- na.omit(USDCOP)

# Gráfico del precio
chartSeries(USDCOP,
            type  = "line",
            theme = chartTheme("white", up.col = "steelblue"),
            TA    = NULL,
            name  = "Tipo de cambio USD/COP (Cierre diario)")

# ─── 2. Cálculo de retornos logarítmicos ──────────────────────────────────── #

# r_t = ln(P_t) - ln(P_{t-1})
USDCOP_R <- diff(log(USDCOP))
USDCOP_R <- na.omit(USDCOP_R)

# Gráfico de retornos
chartSeries(USDCOP_R,
            type  = "line",
            theme = chartTheme("white", up.col = "darkgreen"),
            TA    = NULL,
            name  = "Retornos diarios USD/COP")

# ─── 3. Estimación ARMA(0,0) para obtener residuales ──────────────────────── #
# La media de los retornos financieros se asume cero en el corto plazo.
# El objetivo es aislar los errores para analizar su varianza.

arma00 <- arima(USDCOP_R,
                order        = c(0, 0, 0),
                include.mean = FALSE,
                method       = "ML")

# ─── 4. Diagnóstico de efectos ARCH ───────────────────────────────────────── #

errores <- residuals(arma00)

# 4.1 Gráfico de residuales al cuadrado
chartSeries(errores^2,
            name  = "Residuales al cuadrado - USD/COP",
            theme = chartTheme("white", up.col = "brown"),
            TA    = NULL)

# 4.2 Test formal de efectos ARCH (H0: no hay efectos ARCH)
arch.test(arma00, output = TRUE)
# p-valor < 0.05 rechaza H0 → hay heterocedasticidad condicional

### Se encuentra presencia de efectos ARCH en los residuales de la serie con la presencia de heteroceasticidad y autocorrelación en rezagos cortos. ###

# 4.3 FACP de errores al cuadrado para seleccionar orden
par(mfrow = c(1, 2))
acf(errores^2,  lag.max = 30, main = "FAC  - Residuales²")
pacf(errores^2, lag.max = 30, main = "FACP - Residuales²")
par(mfrow = c(1, 1))

### A partir de las graficas de FAC y FACP se puede identificar la presencia de rezagos significativos, pero se descarta la idea de que sea un modelo ARCH(6). Porque no sostiene su estructura dentro del modelo, por lo que puede ser un efecto estadistico de falso positivo. 
### Se decide estimar el modelo GARCH (1,1) obedeciendo a lo modelos estandar en finanzas y tratando de capturar los efectos ARCH.

# ─── 5. Estimación GARCH(1,1) ─────────────────────────────────────────────── #

spec_garch <- ugarchspec(
  mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  distribution.model = "norm"
)

fit_garch <- ugarchfit(spec = spec_garch, data = USDCOP_R)

# Resultados completos
fit_garch

# Coeficientes
coef(fit_garch)
show(fit_garch)

# Nota
### Al correrse las pruebas de supuestos, se puede encontrar que existen problemas en las pruebas de supuestos.
## 1. Existe una correlación serial en la media, al verse esto en los residuales estandarizados.
## 2. Los parametros no son estables a lo largo de la muestra, presencia de multiples cambios estructurales.
## 3. El efecto conjunto de los shock no son simetricos, cambios positivios o negativos no tienen el mismo efecto en la volatilidad.
## 4. La estructura de distribución de errores en normalidad es rechazada, probar otras distribuciones puede ser recomendable.

# Significancia de los estimadores

# Verificar que existen y tienen la forma correcta
print(fit_garch@fit$matcoef)
print(fit_garch@fit$robust.matcoef)

# Extraer columnas directamente
estimadores <- fit_garch@fit$matcoef[, 1]
SE_std      <- fit_garch@fit$matcoef[, 2]
t_std       <- fit_garch@fit$matcoef[, 3]
p_std       <- fit_garch@fit$matcoef[, 4]

SE_rob      <- fit_garch@fit$robust.matcoef[, 2]
t_rob       <- fit_garch@fit$robust.matcoef[, 3]
p_rob       <- fit_garch@fit$robust.matcoef[, 4]

# ── Función asteriscos─
stars <- function(p) {
  ifelse(p < 0.01,  "***",
         ifelse(p < 0.05,  "**",
                ifelse(p < 0.10,  "*", "")))
}

# Tabla final
tabla_final <- data.frame(
  Estimador = round(estimadores, 8),
  SE_std    = round(SE_std,      8),
  t_std     = round(t_std,       4),
  p_std     = round(p_std,       6),
  sig_std   = stars(p_std),
  SE_rob    = round(SE_rob,      8),
  t_rob     = round(t_rob,       4),
  p_rob     = round(p_rob,       6),
  sig_rob   = stars(p_rob)
)

print(tabla_final)

# ─── 6. Diagnóstico post-estimación ───────────────────────────────────────── #

errores_garch2 <- residuals(fit_garch)^2

par(mfrow = c(1, 2))
acf(errores_garch2,  lag.max = 30, main = "FAC  Residuales² GARCH(1,1)")
pacf(errores_garch2, lag.max = 30, main = "FACP Residuales² GARCH(1,1)")
par(mfrow = c(1, 1))

### Los resultados de las FACP se muestran extremadamente similares entre pre y post GARCH, por lo tanto es posible inferir que la heteroceasticidad no fue del todo capturada.
### El modelo GARCH no captura totalmente la heteroceasticidad de COP/USD
### A pesar de eso, es notable que en gran medida el modelo es funcional pero no perfecto, el rezago no esta diferencialmente alejado de las bandas de confianza.

# ─── 7. Análisis gráfico de la varianza estimada ──────────────────────────── #

var_estimada <- ts(fit_garch@fit$var, start = c(2025, 7, 1), freq = 252)

autoplot(var_estimada) +
  geom_line(color = "royalblue") +
  labs(
    title    = "Varianza condicional estimada - USD/COP",
    subtitle = "Modelo GARCH(1,1)",
    y        = "σ²_t",
    x        = ""
  ) +
  theme_minimal()

# ─── 8. Pronóstico 8 días adelante ─────────────────────────────────────────  #

###──── 8.1 Pronostico de la volatilidad ──────────────────────────────────###
pronostico <- ugarchforecast(fit_garch, n.ahead = 8)
pronostico

# Gráfico del pronóstico de volatilidad (σ_t, no varianza)
plot(pronostico, which = 3)

# Extraer los valores explícitamente
sigma_forecast <- as.data.frame(sigma(pronostico))
print(sigma_forecast)

###─── 8.2 Pronostico con valores de Peso/USD ─────────────────────────────###

# Serie histórica en dataframe 
historico_df <- data.frame(
  Fecha  = index(USDCOP),
  Precio = as.numeric(USDCOP)
)

# Construir el pronóstico desde el último punto histórico 
ultima_fecha  <- max(historico_df$Fecha)
ultimo_precio <- tail(historico_df$Precio, 1)

fechas_forecast <- seq.Date(ultima_fecha, by = "day", length.out = 9)  # día 0 al 8

sigma_h <- as.numeric(sigma(pronostico))   # vector de 8 sigmas

forecast_df <- data.frame(
  Fecha       = fechas_forecast,
  Centro      = ultimo_precio,
  Inferior_95 = c(ultimo_precio, 
                  ultimo_precio * exp(-1.96 * sigma_h * sqrt(1:8))),
  Superior_95 = c(ultimo_precio, 
                  ultimo_precio * exp( 1.96 * sigma_h * sqrt(1:8)))
)

# Gráfico combinado 
ggplot() +
  # Serie histórica
  geom_line(data = historico_df,
            aes(x = Fecha, y = Precio),
            color = "steelblue", linewidth = 0.6) +
  # Banda de predicción
  geom_ribbon(data = forecast_df,
              aes(x = Fecha, ymin = Inferior_95, ymax = Superior_95),
              fill = "tomato", alpha = 0.25) +
  # Línea central del pronóstico
  geom_line(data = forecast_df,
            aes(x = Fecha, y = Centro),
            color = "tomato", linetype = "dashed", linewidth = 0.7) +
  # Punto de unión
  geom_vline(xintercept = ultima_fecha,
             linetype = "dotted", color = "gray40") +
  labs(
    title    = "USD/COP — Serie histórica y pronóstico de volatilidad",
    subtitle = "Intervalo al 95% construido con σ_t del GARCH(1,1) | 8 días adelante",
    x        = "",
    y        = "COP por USD"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(color = "gray40"))

