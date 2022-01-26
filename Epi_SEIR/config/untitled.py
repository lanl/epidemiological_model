fit_disease = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict =  self.fit_out.params)
fit_disease.run_model(disease_name)


model_df = pd.DataFrame(dict(zip(list(fit_disease.initial_states.keys()), fit_disease.model_output.T)))

week_out = model_df.iloc[::7, :].reset_index()
week_out['Dh'] = week_out['Ch'].diff().fillna(0)

data = pd.read_csv('fit_data/rdj_2010_2011_season_peak.csv')

week_out['Time'] = fit_disease.t_eval[::7]

plt.plot(week_out['Time'], data['Dh'], 'ro', label=f'data')
plt.plot(week_out['Time'], week_out['Dh'], 'b-', label= f'fit')