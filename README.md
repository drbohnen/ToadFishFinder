# ToadFishFinder
MATLAB software to identify oyster toadfish calls in passive acoustic data using deep learning. 

Bohnenstiehl, D.R. submitted, Automated Cataloging of Oyster Toadfish (Opsanus tau) Boatwhistle Calls Using Template Matching and Machine Learning   

toadfish_classifier_v4.mat = classifier from the above paper, trained on >20,000 labeled signals from data in western Pamlico Sound, NC.  
This is the recommended classifier. 


Instructions and a brief tutorial are provided on the PerchPicker Wiki! 
https://github.com/drbohnen/ToadFishFinder/wiki


**Overview:** Oyster toadfish (Opsanus tau) represent an ecologically significant species found throughout estuaries along the eastern coast of the United States. While these crevice-dwelling fish can be challenging to observe in their habitats, it is possible to infer their distribution and aspects of their behavior by recording the sounds they produce.  The task of cataloging the distinctive advertisement boatwhistles sounds produced by male toadfish throughout the spring and summer is automated using a multi-step process. Candidate boatwhistles are first identified by template matching using a suite of synthetic spectrogram kernels formed to mimic the two lowest frequency harmonic tones within the boatwhistle. Candidate boatwhistle calls are identified based on the correlation between these kernels and a low-frequency spectrogram of the data. Next, frequency-reassigned spectrogram images of these candidates are formed and input into the pre-trained ResNet-50 convolutional neural network. Finally, the activations from a deep fully connected layer are extracted and passed to a one-vs-all support-vector-machine classifier, which separates boatwhistles from the larger set of candidate signals.  This classifier model was trained and evaluated using a labeled dataset of over 20,000 signals generated over diverse acoustic conditions within Pamlico Sound, North Carolina, USA.  The accompanying software provides an effective and efficient tool to monitor boatwhistles calls, which may facilitate a deeper understanding of the spatial distribution, behavioral patterns, and ecological roles played by oyster toadfish.

