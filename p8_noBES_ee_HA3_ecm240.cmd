Random:setSeed = on
Main:numberOfEvents = 100000         ! number of events to generate
Main:timesAllowErrors = 5          ! how many aborts before run stops

! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showChangedParticleData = off ! list changed particle data
Next:numberCount = 100             ! print message every n events
Next:numberShowInfo = 1            ! print event information n times
Next:numberShowProcess = 1         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

Beams:idA = 11                   ! first beam, e+ = 11
Beams:idB = -11                   ! second beam, e- = -11


! Vertex smearing :
Beams:allowVertexSpread = on
Beams:sigmaVertexX = 9.70e-3   !  13.7 mum / sqrt2
Beams:sigmaVertexY = 25.5E-6   !  36.1 nm / sqrt2
Beams:sigmaVertexZ = 0.64      !  0.64 mm


! 3) Hard process : ZH at 240 GeV
Beams:eCM = 240.  ! CM energy of collision
! HiggsSM:ffbar2HZ = on
Higgs:useBSM = on
HiggsBSM:all = off !Was on
HiggsBSM:ffbar2A3H1 = on


! 4) Settings for the event generation process in the Pythia8 library.
PartonLevel:ISR = on               ! initial-state radiation
PartonLevel:FSR = on               ! final-state radiation

! 5) Non-standard settings; exemplifies tuning possibilities.
25:m0        = 125.0               ! Higgs mass
36:m0        = 5                   ! A3 mass
36:onMode    = off                      ! Switch off A3 decays
36:onIfAny   = 13                            ! Switch on A3 decay to muons
36:mMin      = 1     !set minimum

