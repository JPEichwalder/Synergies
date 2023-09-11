function c3devents = setEventTimes(c3dfile, evtime, label, context, varargin);
  
  ow = 0;%varargin- overwrite; 0 = no = add new events to existing events; 1 = yes = clear existing events
  it = 0;%varargin- intime; 0 = doesn't matter if there are times which are outside of the region; 1= only rite times that are in the region
  while (numel(varargin)>0) %changes values from default to specified input values: *'inputname', value*
    switch (varargin{end-1})
      case 'overwrite'
        ow = varargin{end};
        varargin(end-1:end) = [];
      case 'intime'
        it = varargin{end};
        varargin(end-1:end) = [];
      otherwise
        error ("check input arguments");
    end%switch
  end%while
  
  if it == 1
    endtimec3d    = c3dfile.header.points.lastFrame/c3dfile.header.points.frameRate;
    starttimec3d  = c3dfile.header.points.firstFrame/c3dfile.header.points.frameRate;
    evtime(evtime<starttimec3d) = 0;
    evtime(evtime>endtimec3d) = 0;
    evtime = nonzeros(evtime)';
  end%if
  
  
  
  if ow == 1
    c3dfile.parameters.EVENT.TIMES.DATA(:,:)      = [];
    c3dfile.parameters.EVENT.USED.DATA            = 0;
    c3dfile.parameters.EVENT.CONTEXTS.DATA        = cell(length(evtime),1);
    c3dfile.parameters.EVENT.ICON_IDS.DATA(:)     = [];
    c3dfile.parameters.EVENT.LABELS.DATA          = cell(length(evtime),1);
    c3dfile.parameters.EVENT.DESCRIPTIONS.DATA    = cell(length(evtime),1);
    c3dfile.parameters.EVENT.SUBJECTS.DATA        = cell(length(evtime),1);
    c3dfile.parameters.EVENT.GENERIC_FLAGS.DATA   = [];
  end%if
  
  exiev = c3dfile.parameters.EVENT.USED.DATA;  %nr of already existing events
  newev = length(evtime);
  allev = exiev + newev;  
  
  switch label
    case 'Foot Strike'
      iconID = 1;
      descri = 'The instant the heel strikes the ground';
    case 'Foot Off'
      iconID = 2;
      descri = 'The instant the toe leaves the ground';
    otherwise
      error('check input argument context - has to be Foot Strike or Foot Off')
  end%switch
  

  
  c3dfile.parameters.EVENT.TIMES.DATA(2,exiev+1:allev)      = evtime;
  c3dfile.parameters.EVENT.USED.DATA                        = allev;
  c3dfile.parameters.EVENT.CONTEXTS.DATA(exiev+1:allev)     = {context};
  c3dfile.parameters.EVENT.ICON_IDS.DATA(exiev+1:allev)     = iconID;
  c3dfile.parameters.EVENT.LABELS.DATA(exiev+1:allev)       = {label};
  c3dfile.parameters.EVENT.DESCRIPTIONS.DATA(exiev+1:allev) = {descri};
  c3dfile.parameters.EVENT.SUBJECTS.DATA(exiev+1:allev)     = {c3dfile.parameters.SUBJECTS.NAMES.DATA{1}};
  c3dfile.parameters.EVENT.GENERIC_FLAGS.DATA                = [];
  
  c3devents = c3dfile;
  
end%function

  
  
  
  
  
  
  
  
  
  
  
  