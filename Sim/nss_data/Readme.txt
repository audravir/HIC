Diaa Noureldin, Neil Shephard and Kevin Sheppard, "Multivariate High-Frequency-Based Volatility (HEAVY) Models", 
Journal of Applied Econometrics, forthcoming.

The high-frequency data used in this paper has been obtained from the Trade and Quote (TAQ) database of the New York Stock Exchange (NYSE), 
which can be accessed through the Wharton Research Data Services (WRDS) at http://wrds-web.wharton.upenn.edu/wrds.

The data were cleaned.

The data are organized in 9 Comma Separated Values (.csv) files. Each file contains daily observations for 2242 days. 
These are the trading days on the NYSE during the period 1/2/2001 to 31/12/2009. The following is a description of 
the data files: 

10_dim_daily_return: daily returns for the 10 DJIA stocks (listed below with the ticker symbol). The first column contains the dates, 
the following 10 columns are the Close-to-Close (C-to-C) daily returns, and the last 10 columns are the
Open-to-Close (O-to-C) daily returns.

10_dim_realized_covar: daily realized variances and covariances for the 10 DJIA stocks. The first column contains the dates.
The remaining 55 columns are in the form of the vech of the 10x10 realized covariance matrix where the
stocks are ordered as follows: Bank of America (BAC), JP Morgan (JPM), International Business Machines (IBM), Microsoft (MSFT), 
Exxon Mobil (XOM), Alcoa (AA), American Express (AXP), Du Pont (DD), General Electric (GE) and Coca Cola (KO).  
Applying an inverse vech operator gives the realized covariance matrix.

spy_bac_daily_return: SPY and BAC Close-to-Close (C-to-C) and Open-to-Close (O-to-C) daily returns.
The first column contains the dates, followed by the C-to-C returns and O-to-C, respectively. 

spy_bac_realized_covariance_1min: SPY and BAC realized covariance matrix using 1 minute returns, without subsampling.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance. 

spy_bac_realized_covariance_5min: SPY and BAC realized covariance matrix using 5 minute returns, with subsampling using 1-minute returns.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance. 

spy_bac_realized_covariance_10min: SPY and BAC realized covariance matrix using 10 minute returns, with subsampling using 1-minute returns.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance. 

spy_bac_realized_covariance_15min: SPY and BAC realized covariance matrix using 15 minute returns, with subsampling using 1-minute returns.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance. 

spy_bac_realized_covariance_30min: SPY and BAC realized covariance matrix using 30 minute returns, with subsampling using 1-minute returns.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance. 

spy_bac_realized_kernel: SPY and BAC realized covariance matrix using the realized kernel of 
Ole E. Barndorff-Nielsen, Peter R. Hansen, Asger Lunde and Neil Shephard, 
"Multivariate realised kernels: Consistent positive semi-definite estimators of the covariation of equity prices 
with noise and non-synchronous trading", Journal of Econometrics, Vol. 162, 2011, No. 2, pp. 149-169.
The first column contains the dates, followed by SPY realized variance, BAC realized variance and SPY-BAC realized covariance.

All files are zipped in the file nss-data.zip. 

This directory also includes a supplementary appendix, NSS_appendix.pdf.


   




