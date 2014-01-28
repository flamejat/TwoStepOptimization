function PlotVolumeEvolutionAndPumpSignals( nMainTanks, nAuxiliarTanks, Volumes, signals, perMinuteRate)

	%Volume Evolution
	nTanks = nMainTanks + nAuxiliarTanks;
	figure;
	for i = 1:nMainTanks
		subplot(nTanks,1,i);
		plot(Volumes(:,i));
		axis([0 1500 -5 5])
		hold on
		xlim = get(gca,'xlim');
		plot(xlim, [4 4],'r')
		hold on
		plot(xlim, [0 0],'r')
		box off;
	end

	for i = nMainTanks + (1:nAuxiliarTanks)
		subplot(nTanks, 1, i)
		plot( Volumes( :, i) );
		axis([0 1500 0 5])
		hold on
		plot(xlim, [2.5 2.5], 'r')
		hold on
		plot(xlim, [0 0],'r')
		box off;
	end

	%Pump Signals
	figure;
	for i = 1:nAuxiliarTanks
		subplot(nAuxiliarTanks, 1, i);
		plot(1:(60*24), perMinuteRate/100,'r')
		box off
		hold on;
		pump  = signals(:, i);
		stairs( pump);
	end