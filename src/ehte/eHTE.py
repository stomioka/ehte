# #############################################################################
# Copyright 2023-2024 Sumitomo Pharma America
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# For inquiry about this package, contact Joshua.siegel@us.sumitomo-pharma.com,
# sam.tomioka@us.sumitomo-pharma.com
# 
# This version is compatible with
# numpy               1.22.4
# pandas              2.2.1
#
# #############################################################################

import numpy as np
import pandas as pd
import pyreadstat
import matplotlib.pyplot as plt
import seaborn as sns
import sys

class eHTE_Estimator:
    '''
    indf: is a dataframe with at least TRT01P, TRT01PN, CHG columns. Placebo arm should be coded as TRT01P='Placebo' and TRT01PN = 1. 
    '''
    
    def __init__(self, indf,  n_perms =10000):
        self.indf = indf
        self.n_perms = n_perms +1 # number of permutations
        required_columns = {'TRT01P', 'TRT01PN', 'CHG'}

        if not required_columns.issubset(self.indf.columns):
            sys.exit("Warning: The DataFrame is missing one or more required columns.")

        # Check if 'TRT01P' contains strings
        if not all(self.indf['TRT01P'].apply(lambda x: isinstance(x, str))):
            sys.exit("Warning: 'TRT01P' column should contain only strings.")

        # Check if 'TRT01PN' is integer
        if not all(self.indf['TRT01PN'].apply(lambda x: isinstance(x, int))):
            sys.exit("Warning: 'TRT01PN' column should contain only integers.")

        # Check if 'TRT01P' contains 'Placebo'
        if 'Placebo' not in self.indf['TRT01P'].values:
            sys.exit("Warning: 'TRT01P' column should contain the value 'Placebo'.")

        # Check if 'TRT01PN' is 1 when 'TRT01P' is 'Placebo'
        placebo_rows = self.indf[self.indf['TRT01P'] == 'Placebo']
        if placebo_rows[placebo_rows['TRT01PN'] == 1].empty:
            sys.exit("Warning: 'TRT01PN' should be 1 when 'TRT01P' is 'Placebo'.")

        # Check if 'CHG' is float or integer
        if not all(self.indf['CHG'].apply(lambda x: isinstance(x, (int, float)))):
            sys.exit("Warning: 'CHG' column should contain only integers or floats.")

        else:
            self.trt_dict = indf[['TRT01P', 'TRT01PN']].drop_duplicates().set_index('TRT01PN')['TRT01P'].to_dict()

    def gen_sim(self, trts, seed, nobs_dict, m_dict, s_1):
        '''
        Simulate data for using their observed mean but placebo standard deviation (s_1)
        trts: a list of treatment code in int. e.g. [2,3]
        seed: seed
        s_1: placebo std in float
        m_dict:  mean of each active treatment arm in dict. e.g {2: 1.0, 3: 1.5}
        nobs_dict:  N of each active treatment arm in dict. e.g. {2: 50, 3: 60}
        '''
        np.random.seed(seed)

        data = [
            {'NPERMS': nperms, 'TRT01PN': i, 'PT': pt, 'M': m_dict[i], 'CHG': np.random.normal(m_dict[i], s_1)}
            for nperms in range(1, self.n_perms)
            for i in trts
            for pt in range(1, nobs_dict[i] + 1)
        ]
        return pd.DataFrame(data)

    def calculate_percentile(self, x):
        ranked_values = (x.rank(method='average') ) / (len(x))
        return (ranked_values * 100).astype(int)

    def sas_percentile(self, x, st, en, pctiles=101):
        '''
        Empirical distribution function with averaging (PCTLDEF=5, default option in SAS)
        https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/procstat/procstat_univariate_details14.htm

        reproduce the same code as
        proc univariate data=&indatap. noprint;
          by nperms trt01pn;
          var chg;
          output out=pct_pcb pctlpts=0 to 100 pctlpre=p;
        run;
        '''
        n = len(x)
        x = np.sort(x)
        percentiles = np.linspace(st, en, pctiles)

        result = []
        T =[]
        for  t in percentiles: #0 to 100
            p = t / 100
            np_value = n * p
            j = int(np.floor(np_value))
            g = np_value - j
            #print('t={} p={} np={} j={} g={}'.format(t, p, np_value, j, g))
            if j==0:
                result.append(np.mean([x[0],x[1]]))
                T.append(t)            
            elif g == 0:
                if j<n:
                    result.append(np.mean([x[j-1], x[j]]))
                    T.append(t)
                if j==n:
                    result.append(np.mean([x[j-1]]))
                    T.append(t)
            else:
                if j<n:
                    result.append(x[j])
                    T.append(t)
        return result,T

    def calculate_sigma48(self, indf_p, indf_t, st, en, pctiles, interval95=False):
        '''
        indf_p: placebo dataframe with NPERMS, TRT01PN, CHG
        indf_t: treatment dataframe with NPERMS, TRT01PN, CHG
        intervals: float, intervals of percentiles
        interval95: boolean, central 95% interval
        '''

        if 'NPERMS' not in indf_p.columns:
            indf_p = indf_p.copy()
            indf_p['NPERMS'] = 0
        if 'NPERMS' not in indf_t.columns:
            indf_t = indf_t.copy()
            indf_t['NPERMS'] = 0

        pct_pcb = pd.DataFrame()
        for grp, data in indf_p.groupby(['NPERMS','TRT01PN']).CHG:
            data = data.values
            data.sort()
            res, pctl = self.sas_percentile(data, st, en, pctiles)
            df = pd.DataFrame({'NPERMS':np.repeat(grp[0],len(res)),'PCB_CHG':res,'percentile':pctl})
            pct_pcb = pd.concat([pct_pcb,df])

        pct_pcb['percentile']    =   pct_pcb['percentile'].astype(int)
        pct_trt = pd.DataFrame()
        for grp, data in indf_t.groupby(['NPERMS','TRT01PN']).CHG:
            data = data.values
            data.sort()
            res, pctl = self.sas_percentile(data, st, en, pctiles)
            df = pd.DataFrame({'NPERMS':np.repeat(grp[0],len(res)),'TRT01PN':np.repeat(grp[1],len(res)),'CHG':res,'percentile':pctl})
            pct_trt = pd.concat([pct_trt,df])
        pct_trt['percentile']    =   pct_trt['percentile'].astype(int)    

        #pct_trt['percentile'] = pct_trt.groupby(['NPERMS', 'TRT01PN'])['CHG'].apply(calculate_percentile)

        ite_df = pd.merge(pct_trt, pct_pcb, on=['NPERMS','percentile'], how='left')
        ite_df['ITE'] = ite_df['CHG'] - ite_df['PCB_CHG']
        if interval95 == True: 
            ite_df = ite_df[(ite_df['percentile'] >= 3) & (ite_df['percentile'] <= 97)]

        ehte_df = ite_df.groupby(['NPERMS','TRT01PN'])['ITE'].std().reset_index(name='sigma')

        return ehte_df, ite_df

    def calculate_sigma(self, indf_p, indf_t, st, en, pctiles, interval95):
        '''
        indf_p: placebo dataframe with NPERMS, TRT01PN, CHG
        indf_t: treatment dataframe with NPERMS, TRT01PN, CHG
        intervals: float, intervals of percentiles
        interval95: boolean, central 95% interval
        '''

        if 'NPERMS' not in indf_p.columns:
            indf_p = indf_p.copy()
            indf_p['NPERMS'] = 0
        if 'NPERMS' not in indf_t.columns:
            indf_t = indf_t.copy()
            indf_t['NPERMS'] = 0

        pct_pcb = pd.DataFrame()
        for grp, data in indf_p.groupby(['NPERMS','TRT01PN']).CHG:
            data = data.values
            data.sort()
            res, pctl = self.sas_percentile(data, st, en, pctiles)
            df = pd.DataFrame({'NPERMS':np.repeat(grp[0],len(res)),'PCB_CHG':res,'percentile':pctl})
            pct_pcb = pd.concat([pct_pcb,df])        
        pct_pcb['percentile']    =   pct_pcb['percentile'].astype(int)

        pct_trt = indf_t.copy()

        #percentile_series = pct_trt.groupby(['NPERMS', 'TRT01PN'])['CHG'].apply(self.calculate_percentile)
        pct_trt['percentile'] = pct_trt.groupby(['NPERMS', 'TRT01PN'])['CHG'].transform(lambda x: self.calculate_percentile(x))

        ite_df = pd.merge(pct_trt, pct_pcb, on=['NPERMS','percentile'], how='left')
        ite_df['ITE'] = ite_df['CHG'] - ite_df['PCB_CHG']
        if interval95 == True: 
            ite_df = ite_df[(ite_df['percentile'] >= 3) & (ite_df['percentile'] <= 97)]

        ehte_df = ite_df.groupby(['NPERMS','TRT01PN'])['ITE'].std().reset_index(name='sigma')

        return ehte_df, ite_df

    def cal_pvalue(self, simhte_df, ehte_df, s_1):

        sigmas = {k:v for k,v in zip(ehte_df['TRT01PN'],ehte_df['sigma'])}
        count_ = simhte_df.rename({'sigma':'sim_sigma'}, axis=1).merge(ehte_df, on=['TRT01PN'], how='left')
        count_['ind'] = count_.apply(lambda x: 1 if x.sim_sigma>=x.sigma else 0, axis=1)
        p_values = count_.groupby('TRT01PN').describe()['ind'][['mean']].to_dict()['mean']
        sigma_stds = {k:v for k,v in zip(ehte_df['TRT01PN'], ehte_df.groupby('TRT01PN')['sigma'].apply(lambda x: x / s_1))}
        #print('sigma: {}, eHTE: {}, p-value: {}'.format(ehte_df.sigma[0], sigma_std,p_value))
        for trt_ in list(sigmas.keys()):
            print('treatment: {}, 𝜎: {:2.5f}, eHTE: {:2.3f}, p-value: {:2.3f}'.format(self.trt_dict[int(trt_)], sigmas[trt_], sigma_stds[trt_],p_values[trt_]))
        return sigmas, p_values, sigma_stds

    def eHTE_p(self):
        '''
        indf: a dataframe with TRT01PN, TRT01P, CHG
        '''   
        seeds = [123,456,789,987,654,321,102,103,104,105,106,107,108,109,110]
        summary_stats = self.indf.groupby('TRT01PN').describe()
        display(summary_stats)
        nTRT =self.indf['TRT01PN'].nunique()

        #placebo
        pbo_df = self.indf[(self.indf['TRT01P'].str.lower()=='placebo')].copy()
        placebo = pbo_df.TRT01PN.unique()[0]
        #active trt
        trt_df = self.indf[(self.indf['TRT01P'].str.lower()!='placebo')].copy()    
        trts = np.array([int(i) for i in trt_df['TRT01PN'].unique()])
        trts.sort()

        #for placebo simulation
        result = pbo_df.groupby('TRT01PN').describe()['CHG'][['count','mean','std']]
        result_dict = result.to_dict()
        nobs_dict = {int(key): int(value) for key, value in result.to_dict()['count'].items()}
        m_dict =    {int(key): value for key, value in result.to_dict()['mean'].items()}
        s_ = result.values[0][2] #std of placebo
        #print(s_)
        simp_df = pd.DataFrame()
        simp_df = pd.concat([simp_df,self.gen_sim([placebo],seeds[1],nobs_dict,m_dict,s_)])

        #for trt
        result = trt_df.groupby('TRT01PN').describe()['CHG'][['count','mean']]
        result_dict = result.to_dict()
        nobs_dict = {int(key): int(value) for key, value in result.to_dict()['count'].items()}
        m_dict =    {int(key): value for key, value in result.to_dict()['mean'].items()}
        simt_df = pd.DataFrame()
        simt_df = pd.concat([simt_df,self.gen_sim(trts,seeds[1],nobs_dict,m_dict,s_)])

        #sigma all records
        ehte_df, ite_df = self.calculate_sigma(pbo_df, trt_df, 0,100,101, True)
        simhte_df, _ = self.calculate_sigma(simp_df, simt_df, 0,100,101, True)
        #sigma 48 pctlines
        ehte_48_df, _ = self.calculate_sigma48(pbo_df, trt_df, 3,97,48, False)
        simhte_48_df, _ = self.calculate_sigma48(simp_df, simt_df, 3,97,48, False) 

        print('3-97% percentiles based on actual')
        sigmas, p_values, sigma_stds = self.cal_pvalue(simhte_df, ehte_df, s_)
        print('Fixed 48 percentiles')
        sigmas48, p_values48, sigma_stds48 =  self.cal_pvalue(simhte_48_df, ehte_48_df, s_)

        res_dict = {}
        res_dict['all data'] = {}
        res_dict['48 percentiles'] = {}
        res_dict['all data']['sigmas'] = sigmas
        res_dict['all data']['sigma_stds'] = sigma_stds
        res_dict['all data']['p_values'] = p_values
        res_dict['48 percentiles']['sigmas'] = sigmas
        res_dict['48 percentiles']['sigma_stds'] = sigma_stds
        res_dict['48 percentiles']['p_values'] = p_values

        return ite_df, res_dict

    def eHTEplot(self, ite_df, var_name):
        '''
        assumes TRT01P contains "Placebo"

        '''
        color_lst = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E','#77AC30',"#4DBEEE","#A2142F"]

        trt_lst = set(self.indf[self.indf.TRT01P.str.lower()!='placebo'].TRT01P)

        nn = self.indf.TRT01P.unique().shape[0]
        t_c = color_lst[1:nn]
        cdict = {k:v for k, v in zip(trt_lst,t_c)}
        cdict.update({'Placebo':'#0072BD'}) #placebo is coded as "Placebo" in TPT01P

        fig, [ax0,ax1,ax2] = plt.subplots(figsize=(10,5), ncols=3, gridspec_kw={'width_ratios': [5,3, 1]})

        sns.histplot(data=self.indf, x='CHG', hue='TRT01P',bins=40,multiple='dodge',
                     kde=True,kde_kws={},stat='count', palette=cdict,
                     ax=ax0)
        sns.move_legend(ax0, "upper right",
             title=None, frameon=False)
        ax0.set(xlabel='{}'.format(var_name),ylabel='Frequency')

        sns.ecdfplot(data=self.indf,x='CHG', hue='TRT01P', ax=ax1, palette=cdict)
        ax1.get_legend().remove()
        ax1.set(xlabel='{}'.format(var_name),ylabel='Cummulative Percentile')
        ax1.set_ylim(0,1)
        
        sns.scatterplot(data=ite_df, x='ITE',y='percentile', hue='TRT01P', palette=cdict, ax=ax2)
        ax2.set(xlabel='ITE {}'.format(var_name),ylabel='rank')
        ax2.axvline(x=0, c='black')
        ax2.get_legend().remove()
        ax2.set_ylim(0,100)
        
        for ax in  [ax0,ax1,ax2]:
            ax.spines[['right', 'top']].set_visible(False)
        ax2.spines[['left']].set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)
        plt.show()