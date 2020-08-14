% Yokogawa MEG Reader toolbox for MATLAB
% Version 1.04.01 01-Nov-2011
%
% Files
%   getYkgwData        - Get measured data
%   getYkgwHdrSystem   - Get header of the system information
%   getYkgwHdrChannel  - Get header of the system information
%   getYkgwHdrAcqCond  - Get header of the system information
%   getYkgwHdrEvent    - Get header of the trigger event
%   getYkgwHdrCoregist - Get header of the system information
%   getYkgwHdrDigitize - Get header of the digitization information
%   getYkgwHdrSubject  - Get header of the subject information
%   getYkgwHdrBookmark - Get header of bookmark
%   getYkgwHdrSource   - Get header of the source information
%   getYkgwMriHdr      - Get header of mri
%   getYkgwVersion     - Get toolbox version
%
%------------------------------------------------------------
%                          History
%------------------------------------------------------------
% R1.04.01 : 2011.11.01
%    - License Agreement was updated.
%
% R1.04.00 : 2011.05.06
%  [ getYkgwHdrAcqCond ]
%    - Add multi-trigger information
%  [ getYkgwHdrEvent ]
%    - support only EvokedRaw or EvokedAve
%  [ getYkgwHdrChannel ]
%    - remove the end of unnecessary spaces for channel name
%  [ getYkgwHdrDigitize ]
%    - fix .meg2digitizer .digitizer2meg
%  [ getYkgwMriHdr ]
%    - add HPI(marker) information
%
% R1.03.01 : 2011.04.13
%  [ getYkgwData ]
%    - Support internal function changes
%
% R1.03.00.01 : 2011.03.25
%  [ all functions ]
%    - Pcodes were generated by MATLAB R2007b.
%
% R1.03.00 : 2011.03.03
%  [ getYkgwData ]
%    - Support MegLaboratory R1.4.8 specific calibration for non-MEG channels
%  [ getYkgwHdrSource ]
%    - Structure fields which had not been used were removed.
%  [ getYkgwHdrBookmark ]
%    - Structure fields which had not been used were removed.
%
% R1.02.00 : 2011.02.14
%  [ all functions ]
%    - 1st argument (specifying a file) was modified from file ID to file path.
%  [ getYkgwData ]
%    - Fix EvokedRaw error if sample_length argument is omitted.
%  [ getYkgwHdrCoregist ]
%    - Add HPI(Head Position Indicator) fields (position and label)
%
% R1.01.00 : 2011.01.24
%  [ getYkgwData ]
%    - Not apply averaged gain to reference channels if system id is 1000 later. 
%
% R1.00.00 : 2010.06.24
%  - First release
% 
%-----------------------------------------------------------------------
% License Agreement							
%
% Copyright (c) 2010-2011 YOKOGAWA Electric Corporation.
% All rights reserved.
% 
% Yokogawa MEG Reader Toolbox is distributed subject to the following license conditions:
% 
% SOFTWARE LICENSE AGREEMENT*
% Software: Yokogawa MEG Reader Toolbox
% 
% 1. The "software", below, refers to Yokogawa MEG Reader Toolbox (in binary form and accompanying documentation). 
% Each licensee is addressed as "you" or "Licensee."
% 
% 2. The copyright holders shown above and their third-party licensors hereby grant Licensee a royalty-free nonexclusive license,
% subject to the limitations stated herein.
% 
% 3. You may make a copy or copies of the software for use within your organization, if you meet the following conditions:
% Copies must include the copyright notice and this software License Agreement in the documentation
% and/or other materials provided with the copy.
% 
% 4. You may not modify a copy or copies of the software or any portion of it. 
% 
% 5. This software and data processed by this software may not be used for clinical and/or diagnostic purposes.
% 
% 6. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND.
% THE COPYRIGHT HOLDERS, THEIR THIRD PARTY LICENSORS, AND THEIR EMPLOYEES:
% (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED 
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE OR NON-INFRINGEMENT,
% (2) DO NOT ASSUME ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
% COMPLETENESS, OR USEFULNESS OF THE SOFTWARE,
% (3) DO NOT REPRESENT THAT USE OF THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, 
% (4) DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION UNINTERRUPTED, 
% THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL BE CORRECTED.
% 
% 7. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT HOLDERS, THEIR THIRD PARTY LICENSORS, 
% OR THEIR EMPLOYEES: BE LIABLE TO YOU AND/OR ANY PARTY FOR ANY DIRECT, INDIRECT, SPECIAL, 
% INCIDENTAL OR OTHER CONSEQUENTIAL DAMAGES, OR FOR EXEMPLARY, SPECIAL, 
% PUNITIVE OR SIMILAR DAMAGES OF ANY KIND, WHETHER BASED ON CONTRACT, 
% STRICT LIABILITY, TORT, WARRANTY (EXPRESS OR IMPLIED), OR ANY OTHER LEGAL GROUNDS 
% FOR ANY USE OF THIS SOFTWARE, INCLUDING WITHOUT LIMITATION, ANY LOST PROFITS, 
% LOSS OF INCOME, BUSINESS INTERRUPTION, LOSS OF USE, LOSS OR DESTRUCTION OF PROGRAMS 
% OR OTHER DATA, LOSS OF AVAILABILITY AND THE LIKE ON YOUR INFORMATION HANDLING SYSTEM, 
% OR OTHERWISE, EVEN IF THE COPYRIGHT HOLDERS HAS BEEN EXPRESSLY ADVISED OF THE POSSIBILITY 
% OF SUCH DAMAGES.
% 
% 8. The copyright holders reserve the right to change, or temporarily or permanently withdraw, 
% information and/or contents contained in its own software, or to suspend or discontinue the services 
% provided through this software at any time, without prior notice. 
% The copyright holders shall not be held liable for damage of any kind sustained by you for any reason,
% as a result of the copyright holders changing, or temporarily or permanently withdrawing, such information.
% 
% *This License Agreement is subject to change without notice.