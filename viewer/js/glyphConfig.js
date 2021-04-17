

function glyphSettings()
{
    var out = [

        {gene: 'Snca',          taxonomy: 'in_general',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Cplx2',         taxonomy: 'in_general',  glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Lhx6',          taxonomy: 'in_general',  glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Col25a1',       taxonomy: 'in_general',  glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Pnoc',          taxonomy: 'in_general',  glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Rab3c',         taxonomy: 'in_general',  glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: 'Gad1',          taxonomy: 'in_general',  glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Slc6a1',        taxonomy: 'in_general',  glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Th',            taxonomy: 'sst',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Crhbp',         taxonomy: 'sst',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Sst',           taxonomy: 'sst',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Npy',           taxonomy: 'sst',         glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Synpr',         taxonomy: 'sst',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Chodl',         taxonomy: 'sst',         glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Cort',          taxonomy: 'sst',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Reln',          taxonomy: 'sst',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Serpini1',      taxonomy: 'sst',         glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: 'Satb1',         taxonomy: 'sst',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Grin3a',        taxonomy: 'sst',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Tac1',          taxonomy: 'pvalb',       glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Pvalb',         taxonomy: 'pvalb',       glyphSymbol: '*',  glyphName: 'asterisk'},
       {gene: 'Kcnip2',        taxonomy: 'pvalb',       glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Thsd7a',        taxonomy: 'pvalb',       glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Cox6a2',        taxonomy: 'pvalb',       glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Chrm2',         taxonomy: 'pvalb',       glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Id2',           taxonomy: 'ngf',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Hapln1',        taxonomy: 'ngf',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Gabrd',         taxonomy: 'ngf',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Cryab',         taxonomy: 'ngf',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Kit',           taxonomy: 'ngf',         glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Ndnf',          taxonomy: 'ngf',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Nos1',          taxonomy: 'ngf',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Lamp5',         taxonomy: 'ngf',         glyphSymbol: '>',  glyphName: 'triangleRight'},
       {gene: 'Cplx3',         taxonomy: 'ngf',         glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Cadps2',        taxonomy: 'cxcl14',      glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Cxcl14',        taxonomy: 'cxcl14',      glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Ntng1',         taxonomy: 'cxcl14',      glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Cpne5',         taxonomy: 'cxcl14',      glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Rgs12',         taxonomy: 'cxcl14',      glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Sncg',          taxonomy: 'cnr1',        glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Cnr1',          taxonomy: 'cnr1',        glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Cck',           taxonomy: 'cnr1',        glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Trp53i11',      taxonomy: 'cnr1',        glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Sema3c',        taxonomy: 'cnr1',        glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Syt6',          taxonomy: 'cnr1',        glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Yjefn3',        taxonomy: 'cnr1',        glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Rgs10',         taxonomy: 'cnr1',        glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Nov',           taxonomy: 'cnr1',        glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: 'Kctd12',        taxonomy: 'cnr1',        glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Slc17a8',       taxonomy: 'cnr1',        glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Tac2',          taxonomy: 'vip',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Npy2r',         taxonomy: 'vip',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Calb2',         taxonomy: 'vip',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Htr3a',         taxonomy: 'vip',         glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Slc5a7',        taxonomy: 'vip',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Penk',          taxonomy: 'vip',         glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Pthlh',         taxonomy: 'vip',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Vip',           taxonomy: 'vip',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Crh',           taxonomy: 'vip',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Qrfpr',         taxonomy: 'vip',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Zcchc12',       taxonomy: 'less_active', glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Calb1',         taxonomy: 'less_active', glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Vsnl1',         taxonomy: 'less_active', glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Tmsb10',        taxonomy: 'less_active', glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Rbp4',          taxonomy: 'less_active', glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Fxyd6',         taxonomy: 'less_active', glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: '6330403K07Rik', taxonomy: 'less_active', glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: 'Scg2',          taxonomy: 'less_active', glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Gap43',         taxonomy: 'less_active', glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Nrsn1',         taxonomy: 'less_active', glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Gda',           taxonomy: 'pc_or_in',    glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Bcl11b',        taxonomy: 'pc_or_in',    glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Rgs4',          taxonomy: 'pc_or_in',    glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Slc24a2',       taxonomy: 'pc_or_in',    glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Lphn2',         taxonomy: 'pc_or_in',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Adgrl2',        taxonomy: 'pc_or_in',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Map2',          taxonomy: 'pc_or_in',    glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Prkca',         taxonomy: 'pc_or_in',    glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Cdh13',         taxonomy: 'pc_or_in',    glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Atp1b1',        taxonomy: 'pc_or_in',    glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Pde1a',         taxonomy: 'pc_or_in',    glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: 'Calm2',         taxonomy: 'pc_or_in',    glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Sema3e',        taxonomy: 'pc_or_in',    glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Nrn1',          taxonomy: 'pc',          glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Pcp4',          taxonomy: 'pc',          glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Rprm',          taxonomy: 'pc',          glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Enpp2',         taxonomy: 'pc',          glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Rorb',          taxonomy: 'pc',          glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Rasgrf2',       taxonomy: 'pc',          glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Wfs1',          taxonomy: 'pc',          glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene: 'Fos',           taxonomy: 'pc',          glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Plcxd2',        taxonomy: 'pc',          glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Crym',          taxonomy: 'pc',          glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene: '3110035E14Rik', taxonomy: 'pc',          glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene: 'Foxp2',         taxonomy: 'pc',          glyphSymbol: 'p',  glyphName: 'star5'},
        {gene: 'Pvrl3',         taxonomy: 'pc',          glyphSymbol: 'h',  glyphName: 'star6'},
        {gene: 'Neurod6',       taxonomy: 'pc2',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Nr4a2',         taxonomy: 'pc2',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Cux2',          taxonomy: 'pc2',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Kcnk2',         taxonomy: 'pc2',         glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Arpp21',        taxonomy: 'pc2',         glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Enc1',          taxonomy: 'pc2',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Fam19a1',       taxonomy: 'pc2',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene: 'Vim',           taxonomy: 'non_neuron',  glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene: 'Slc1a2',        taxonomy: 'non_neuron',  glyphSymbol: '.',  glyphName: 'point'},
        {gene: 'Pax6',          taxonomy: 'non_neuron',  glyphSymbol: 's',  glyphName: 'square'},
        {gene: 'Plp1',          taxonomy: 'non_neuron',  glyphSymbol: 'x',  glyphName: 'cross'},
        {gene: 'Mal',           taxonomy: 'non_neuron',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene: 'Aldoc',         taxonomy: 'non_neuron',  glyphSymbol: 'o',  glyphName: 'circle'},
        {gene: 'Actb',          taxonomy: 'non_neuron',  glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene: 'Sulf2',         taxonomy: 'non_neuron',  glyphSymbol: 'p',  glyphName: 'star5'},

        {gene: 'Generic',       taxonomy: 'generic',     glyphSymbol: 'o',  glyphName: 'circle'},

        ];

    return out
}

//create color ramp.
function glyphColor(y) {
    return y === 'non_neuron' ? '#FFFFFF' : //hsv: [0 0 1]);
        y === 'pc_or_in' ? '#407F59' :      //hsv: [.4 .5 .5]);
            y === 'less_active' ? '#96B38F' :   //hsv: [.3 .2 .7]);
                y === 'pc' ? '#00FF00' :            //hsv: [1/3 1 1]);
                    y === 'pc2' ? '#44B300' :           //hsv: [.27 1 .7]);
                        y === 'in_general' ? '#0000FF' :    //hsv: [2/3 1 1]);
                            y === 'sst' ? '#00B3FF' :           //hsv: [.55 1 1]);
                                y === 'pvalb' ? '#5C33FF' :         //hsv: [.7 .8 1]);
                                    y === 'ngf' ? '#FF00E6' :           //hsv: [.85 1 1]);
                                        y === 'cnr1' ? '#FF0000' :          //hsv: [ 1 1 1]);
                                            y === 'vip' ? '#FFC700' :           //hsv: [.13 1 1]);
                                                y === 'cxcl14' ? '#995C00' :        //hsv: [.1 1 .6]);
                                                    '#FD6A02';
}