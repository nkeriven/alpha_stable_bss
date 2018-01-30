% test script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The  must install the following toolboxes:
%   The bss_eval toolbox (http://bass-db.gforge.inria.fr/bss_eval/)
%   The sketchml toolbox (http://sketchml.gforge.inria.fr)
% All toolboxes must be in matlab path.
%
% The data used in the paper are adapted from the QUASI dataset and TIMIT
% dataset.
% The TIMIT dataset is provided
% the QUASI dataset can be downloaded at http://www.tsi.telecom-paristech.fr/aao/en/2012/03/12/quasi/

clear all
close all
rng(0)

%%%% Datasets
% QUASI DATASET
QUASI_PATH = 'C:\Users\nkeriven\Documents\local_data\QUASI_QISS\';
QUASI_NAMES = { % list of songs in the QUASI dataset
    '01_alexq-carol_of_the_bells__snip_0-30__'
    '02_another_dreamer-one__snip_69-99__'
    '03_carl_leth-the_world_is_under_attack__snip_44-74__'
    '05_fort_minor-remember_the_name__snip_54-84__'
    '06_glen_philips-the_spirit_of_shackleton__snip_163-193__'
    '07_jims_big_ego-mix_tape__snip_25-55__'
    '08_nine_inch_nails-good_soldier__snip_104-134__'
    '09_shannon_hurley-sunrise__snip_62-92__'
    '10_ultimate_nz_tour__snip_43-73__'
    '11_vieux_farka-ana__snip_120-150__'
    };
QUASI_INSTRUS = { % ordered list of instruments for each song
    {'rythm_gtr','bass','drums','lead_gtr'} % 01
    {'drums','git','bass','lead_vocal'} % 02
    {'drums','piano_loop','piano','elec_gtr','speech','synth'} % 03 (lot of silence)
    {'bass','voc','drums','strings','samples'} % 05
    {'drums','lead_voca','organ','bass','ac_gtr1','pads'} % 06
    {'drums','bass','gtr1','vocs'} % 07
    {'bass','drums','lead_vocal','gtr','vibes','keyboards'} % 08
    {'acoustik_gtrs','cello','drums','bass','voc','piano','elec_gtrs'} % 09
    {'synth','drums','gtr','vocs','bass'} % 10
    {'clav','bass','drums','git','organ','vox','wind'} % 11
    };

% TIMIT DATASET
TIMIT_PATH = 'C:\Users\nkeriven\Documents\local_data\long_TIMIT\';
timit_perm=randperm(10); % TIMIT tests use k tracks randomly selected


%%%% Parameters of the experiment
data_to_use = 'quasi'; % 'quasi' or 'timit'
song_id = 1; % for QUASI: choose song
start_time = 0; % (s) Start time when loading sounds
end_time = 30; % (s) End time when loading sounds 

fs = 16000; % (Hz) Frequency of sampling
tw = 32; % (ms) STFT analysis window
ov = 50; % STFT overlap (%)

M = 2; % Number of microphones
K = 3; % Number of sources

type_estimator = 'alpha'; % 'alpha', 'em_gmm', 'sawada', 'sk_gmm'
sk_size = 200; % sketch size for sketching methods
%%%% the sketching estimators can be long (especially alpha)

% Mixing
max_delay = 20;
max_gain = 5;
min_delta = 1;


% ========================= load data
source_path = cell(K,1);
for k = 1:K
    switch data_to_use
        case 'quasi'
            source_path{k}=[QUASI_PATH,QUASI_NAMES{song_id},QUASI_INSTRUS{song_id}{k},'.wav'];
        case 'timit'
            source_path{k}=[TIMIT_PATH,int2str(timit_perm(k)),'.wav'];
        otherwise
            error('Unrecognized data.')
    end
end

%%%% Load sources
source = squeeze(load_sources(source_path,fs,start_time,end_time,[])); %L*K
L = size(source,1);

%%%% Mix sources with random delays and gains
% mix (L*M), s_gt (L*M*K), delay (M*K), gain (M*K)
[mix,s_gt,delay,gain] = random_mix(source,M,...
                                   max_delay,max_gain,min_delta);

%%%% Short time Fourier domain signals
S_img = mult_stft(s_gt,fs,tw,ov); % F*T*M*K
A_gt = get_mixing_matrices(S_img); % F*M*K : True mixing matrices
X = mult_stft(mix,fs,tw,ov); % F*T*M
[F,T,~] = size(X);

%%%% Estimate masks
[mask,A,noise,alpha] = estimate_masks(X,K,'type_estimator',type_estimator, 'sk_size', sk_size, 'min_alpha', 1); % F*T*K
[mask_perm,A_perm] = oracle_permute(mask,X,S_img,A,alpha,noise); % F*T*K

%%%% Oracle masks to compare
mask_oracle = oracle_masks(S_img); % F*T*K

%%%% Apply masks
S_estim = bsxfun(@times,X,reshape(mask_perm,[F,T,1,K])); % F*T*M*K
S_oracle = bsxfun(@times,X,reshape(mask_oracle,[F,T,1,K])); % F*T*M*K

%%%% Reconstruct time domain signals
s_estim = mult_istft(S_estim,fs,tw,ov,L); % L*M*K
s_oracle = mult_istft(S_oracle,fs,tw,ov,L); % L*M*K
s_mix = repmat(mix,[1,1,K]); % L*M*K

%%%% Source separation evaluation (use bss_eval toolbox)
fprintf(1,'==== Separation results ====\n');
% Oracle
[SDR_oracle,~,SIR_oracle]=bss_eval_images(permute(s_oracle,[3,1,2]),...
                                           permute(s_gt,[3,1,2]));
fprintf(1,'  Oracle : SDR : %.2f+-%.2f  SIR : %.2f+-%.2f\n',...
           mean(SDR_oracle),std(SDR_oracle),mean(SIR_oracle),std(SIR_oracle));
% Estimated
[SDR_estim,~,SIR_estim]=bss_eval_images(permute(s_estim,[3,1,2]),...
    permute(s_gt,[3,1,2]));
fprintf(1,'  Estimated : SDR : %.2f+-%.2f  SIR : %.2f+-%.2f\n',...
           mean(SDR_estim),std(SDR_estim),mean(SIR_estim),std(SIR_estim));
% Mixture (baseline)
[SDR_mix,~,SIR_mix]=bss_eval_images(permute(s_mix,[3,1,2]),...
                                     permute(s_gt,[3,1,2]));
fprintf(1,'  Mixture : SDR : %.2f+-%.2f  SIR : %.2f+-%.2f\n',...
           mean(SDR_mix),std(SDR_mix),mean(SIR_mix),std(SIR_mix));

%%%% compute MER: estimation matrix error
[MER] = bss_eval_mix(permute(A_perm,[2 3 1]), permute(A_gt,[2 3 1]));
fprintf(1,'  MER : %.2f+-%.2f\n',nanmean(MER),nanstd(MER));


