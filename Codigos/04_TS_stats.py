"""
Clase para análisis de series temporales con cálculo de AUC
Autor: Versión refactorizada
"""

import numpy as np
import pandas as pd
from scipy.signal import welch
from scipy import integrate
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.stattools import acf
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field


@dataclass
class AnalysisConfig:
    """Configuración para el análisis de series temporales"""
    fs: float = 1.0  # 1 observación por semana
    nperseg: int = 52  # ventana anual
    freq_muestreo: int = 52  # frecuencia de muestreo anual (semanas)
    acf_lags: List[int] = field(default_factory=lambda: [26, 52, 78, 104, 156])
    qtest_lags: List[int] = field(default_factory=lambda: [1, 2, 3, 26, 52, 104, 156])
    season_map: Dict[int, str] = field(default_factory=lambda: {
        1: 'Winter', 2: 'Spring', 3: 'Summer', 4: 'Fall'
    })


class TimeSeriesAnalyzer:
    """
    Analizador de series temporales para datos de coherencia SAR
    
    Calcula:
    - Área bajo la curva (AUC) del año medio por banda
    - Estadísticas por estación
    - Función de autocorrelación (ACF)
    - Test de Ljung-Box (Q-test)
    - Índices del periodograma
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        """
        Inicializar el analizador
        
        Args:
            config: Configuración personalizada. Si es None, usa valores por defecto
        """
        self.config = config or AnalysisConfig()
        self.data: Optional[pd.DataFrame] = None
        self.results: List[Dict] = []
        
    def load_data(self, input_file: str) -> pd.DataFrame:
        """
        Cargar datos desde archivo CSV
        
        Args:
            input_file: Ruta al archivo CSV
            
        Returns:
            DataFrame con los datos cargados
            
        Raises:
            FileNotFoundError: Si el archivo no existe
        """
        if not Path(input_file).exists():
            raise FileNotFoundError(f"El archivo {input_file} no existe")
        
        self.data = pd.read_csv(input_file)
        self.data['fecha'] = pd.to_datetime(self.data['fecha'])
        
        print(f"✅ Datos cargados desde: {input_file}")
        print(f"   Registros: {len(self.data)}")
        print(f"   Columnas: {self.data.columns.tolist()}")
        print(f"   Bandas: {self.data['band'].unique()}")
        
        return self.data
    
    def _calculate_auc_yearly_mean(self, data_band: pd.DataFrame) -> Dict[str, float]:
        """
        Calcular AUC del año medio
        
        Args:
            data_band: DataFrame filtrado por banda
            
        Returns:
            Dict con métricas de AUC
        """
        # Obtener semana del año
        data_band = data_band.copy()
        data_band['week_of_year'] = data_band['fecha'].dt.isocalendar().week
        
        # Agrupar por semana del año
        yearly_mean = (
            data_band
            .groupby('week_of_year')['coherencia']
            .mean()
            .sort_index()
        )
        
        if len(yearly_mean) == 0:
            return {
                'auc_yearly_mean': np.nan,
                'auc_yearly_mean_normalized': np.nan,
                'mean_yearly_value': np.nan,
                'min_yearly_value': np.nan,
                'max_yearly_value': np.nan,
                'std_yearly_value': np.nan
            }
        
        # Calcular AUC
        x = np.arange(len(yearly_mean))
        y = yearly_mean.values
        auc = integrate.trapezoid(y, x)
        
        return {
            'auc_yearly_mean': auc,
            'auc_yearly_mean_normalized': auc / len(yearly_mean),
            'mean_yearly_value': np.mean(y),
            'min_yearly_value': np.min(y),
            'max_yearly_value': np.max(y),
            'std_yearly_value': np.std(y)
        }
    
    def _calculate_seasonal_stats(self, data_band: pd.DataFrame) -> Dict[str, float]:
        """
        Calcular estadísticas por estación
        
        Args:
            data_band: DataFrame filtrado por banda
            
        Returns:
            Dict con media y desviación estándar por estación
        """
        data_band = data_band.copy()
        data_band['season'] = data_band['fecha'].dt.quarter.map(self.config.season_map)
        
        stats = {}
        for season in ['Winter', 'Spring', 'Summer', 'Fall']:
            season_data = data_band[data_band['season'] == season]['coherencia']
            stats[f'mean_{season}'] = season_data.mean() if len(season_data) > 0 else np.nan
            stats[f'std_{season}'] = season_data.std() if len(season_data) > 0 else np.nan
        
        return stats
    
    def _calculate_acf(self, ts: pd.Series) -> Dict[str, float]:
        """
        Calcular función de autocorrelación en lags específicos
        
        Args:
            ts: Serie temporal
            
        Returns:
            Dict con valores ACF por lag
        """
        try:
            acf_values = acf(ts.values, nlags=max(self.config.acf_lags), fft=True)
            
            return {
                f'acf_lag_{lag}': acf_values[lag] if lag < len(acf_values) else np.nan
                for lag in self.config.acf_lags
            }
        except Exception as e:
            print(f"⚠️  Error calculando ACF: {e}")
            return {f'acf_lag_{lag}': np.nan for lag in self.config.acf_lags}
    
    def _calculate_ljung_box(self, ts: pd.Series) -> Dict[str, float]:
        """
        Calcular test de Ljung-Box (Q-test)
        
        Args:
            ts: Serie temporal
            
        Returns:
            Dict con estadísticos y p-valores
        """
        try:
            lb_result = acorr_ljungbox(ts.values, lags=self.config.qtest_lags, return_df=True)
            
            stats = {}
            for lag in self.config.qtest_lags:
                if lag in lb_result.index:
                    stats[f'qtest_stat_{lag}'] = lb_result.loc[lag, 'lb_stat']
                    stats[f'qtest_pvalue_{lag}'] = lb_result.loc[lag, 'lb_pvalue']
                else:
                    stats[f'qtest_stat_{lag}'] = np.nan
                    stats[f'qtest_pvalue_{lag}'] = np.nan
            
            return stats
        except Exception as e:
            print(f"⚠️  Error calculando Ljung-Box: {e}")
            return {
                f'qtest_stat_{lag}': np.nan for lag in self.config.qtest_lags
            } | {
                f'qtest_pvalue_{lag}': np.nan for lag in self.config.qtest_lags
            }
    
    def _calculate_periodogram_indices(self, ts: pd.Series) -> Dict[str, float]:
        """
        Calcular índices del periodograma
        
        Args:
            ts: Serie temporal
            
        Returns:
            Dict con índices espectrales
        """
        try:
            # Calcular periodograma con Welch
            freqs, power = welch(
                ts.values,
                fs=self.config.fs,
                detrend='constant',
                scaling='density'
            )
            
            # Convertir frecuencias a periodos
            mask = freqs > 0
            periodos = 1.0 / freqs[mask]
            power = power[mask]
            
            # Ciclos de interés
            ciclos = np.array([
                self.config.freq_muestreo / 3,
                self.config.freq_muestreo / 2,
                self.config.freq_muestreo
            ])
            
            # Encontrar posiciones de los ciclos
            pos_ciclos = [np.argmin(np.abs(periodos - ciclo)) for ciclo in ciclos]
            bands_power = power[pos_ciclos]
            
            # Calcular índices
            max_power = np.max(bands_power)
            mean_power = np.mean(bands_power)
            max_idx = np.argmax(bands_power)
            
            power_ciclo_anual_adelante = np.sum(power[pos_ciclos[2]:])
            power_plurianual = np.sum(power[:pos_ciclos[2]])
            power_total = np.sum(power)
            
            return {
                'seasonality_mode': ciclos[max_idx],
                'fisher_kappa': max_power / mean_power if mean_power > 0 else np.nan,
                'seasonality_stability': (
                    max_power / power_ciclo_anual_adelante 
                    if power_ciclo_anual_adelante > 0 else np.nan
                ),
                'plurianual_cycles': (
                    power_plurianual / power_total 
                    if power_total > 0 else np.nan
                ),
                'seasonality_amplitude': max_power
            }
            
        except Exception as e:
            print(f"⚠️  Error calculando periodograma: {e}")
            return {
                'seasonality_mode': np.nan,
                'fisher_kappa': np.nan,
                'seasonality_stability': np.nan,
                'plurianual_cycles': np.nan,
                'seasonality_amplitude': np.nan
            }
    
    def analyze_parcel_band(
        self, 
        fid: int, 
        band: str, 
        sp_ppal: str
    ) -> Optional[Dict]:
        """
        Analizar una parcela específica para una banda
        
        Args:
            fid: ID de la parcela
            band: Banda a analizar
            sp_ppal: Especie principal
            
        Returns:
            Dict con todas las métricas calculadas, o None si no hay datos suficientes
        """
        if self.data is None:
            raise ValueError("Primero debes cargar los datos con load_data()")
        
        # Filtrar datos
        data_fid_band = self.data[
            (self.data['PARCELA'] == fid) &
            (self.data['band'] == band)
        ].copy()
        
        if len(data_fid_band) == 0:
            return None
        
        # Serie temporal completa
        ts = (
            data_fid_band
            .groupby('fecha')['coherencia']
            .mean()
            .dropna()
            .sort_index()
        )
        
        if len(ts) < self.config.nperseg:
            return None
        
        # Inicializar resultado
        result = {
            'PARCELA': fid,
            'SP_PPAL': sp_ppal,
            'band': band
        }
        
        # Calcular todas las métricas
        result.update(self._calculate_auc_yearly_mean(data_fid_band))
        result.update(self._calculate_seasonal_stats(data_fid_band))
        result.update(self._calculate_acf(ts))
        result.update(self._calculate_ljung_box(ts))
        result.update(self._calculate_periodogram_indices(ts))
        
        return result
    
    def analyze_all(self, verbose: bool = True) -> pd.DataFrame:
        """
        Analizar todas las parcelas y bandas
        
        Args:
            verbose: Mostrar progreso durante el análisis
            
        Returns:
            DataFrame con todos los resultados
        """
        if self.data is None:
            raise ValueError("Primero debes cargar los datos con load_data()")
        
        self.results = []
        bands = self.data['band'].unique()
        parcelas = self.data['PARCELA'].unique()
        total_parcelas = len(parcelas)
        
        print("=" * 80)
        print("INICIANDO ANÁLISIS DE SERIES TEMPORALES")
        print("=" * 80)
        
        for idx, fid in enumerate(parcelas, 1):
            if verbose and idx % 50 == 0:
                print(f"Procesando parcela {idx}/{total_parcelas}...")
            
            sp_ppal = self.data[self.data['PARCELA'] == fid]['SP_PPAL'].iloc[0]
            
            for band in bands:
                result = self.analyze_parcel_band(fid, band, sp_ppal)
                if result is not None:
                    self.results.append(result)
        
        return self.get_results()
    
    def get_results(self) -> pd.DataFrame:
        """
        Obtener resultados como DataFrame
        
        Returns:
            DataFrame con todos los resultados
        """
        if not self.results:
            return pd.DataFrame()
        
        df = pd.DataFrame(self.results)
        
        # Reordenar columnas
        cols = ['PARCELA', 'SP_PPAL', 'band'] + [
            col for col in df.columns 
            if col not in ['PARCELA', 'SP_PPAL', 'band']
        ]
        
        return df[cols]
    
    def save_results(self, output_file: str) -> None:
        """
        Guardar resultados en CSV
        
        Args:
            output_file: Ruta del archivo de salida
        """
        df = self.get_results()
        df.to_csv(output_file, index=False)
        
        print("\n" + "=" * 80)
        print("✅ ANÁLISIS COMPLETADO")
        print("=" * 80)
        print(f"Archivo guardado en: {output_file}")
        print(f"Total de registros: {len(df)}")
        print(f"\nColumnas generadas ({len(df.columns)}):")
        for i, col in enumerate(df.columns, 1):
            print(f"  {i:2d}. {col}")
    
    def summary_by_band(self) -> pd.DataFrame:
        """
        Generar resumen de AUC por banda
        
        Returns:
            DataFrame con estadísticas descriptivas por banda
        """
        df = self.get_results()
        if df.empty:
            return pd.DataFrame()
        
        auc_cols = ['auc_yearly_mean', 'auc_yearly_mean_normalized', 'mean_yearly_value']
        return df.groupby('band')[auc_cols].describe()


# ============================================================================
# FUNCIONES DE UTILIDAD
# ============================================================================

def run_standard_analysis(
    input_file: str,
    output_file: str,
    config: Optional[AnalysisConfig] = None
) -> Tuple[TimeSeriesAnalyzer, pd.DataFrame]:
    """
    Ejecutar análisis estándar completo
    
    Args:
        input_file: Ruta del archivo CSV de entrada
        output_file: Ruta del archivo CSV de salida
        config: Configuración personalizada (opcional)
        
    Returns:
        Tupla con (analyzer, results_df)
    """
    analyzer = TimeSeriesAnalyzer(config)
    analyzer.load_data(input_file)
    results = analyzer.analyze_all()
    analyzer.save_results(output_file)
    
    # Mostrar resumen
    print(f"\n📊 Resumen de AUC por banda:")
    summary = analyzer.summary_by_band()
    print(summary)
    
    print(f"\nPrimeras 5 filas:")
    print(results.head())
    
    return analyzer, results


# ============================================================================
# EJEMPLO DE USO
# ============================================================================

if __name__ == "__main__":
    # Configuración
    input_file = r"E:\SAR_UVa\CSV'S\valores_coherencia_extremadura_modificado.csv"
    output_file = 'estadisticas_series_temporales_.csv'
    
    # Opción 1: Usar función de utilidad (más simple)
    analyzer, results = run_standard_analysis(input_file, output_file)
    
    # Opción 2: Uso personalizado (más control)
    # config = AnalysisConfig(
    #     acf_lags=[26, 52, 104],
    #     qtest_lags=[1, 26, 52]
    # )
    # analyzer = TimeSeriesAnalyzer(config)
    # analyzer.load_data(input_file)
    # results = analyzer.analyze_all()
    # analyzer.save_results(output_file)
    
    # Análisis de una parcela específica
    # result_single = analyzer.analyze_parcel_band(fid=1, band='VV', sp_ppal='Pinus')
    # print(result_single)