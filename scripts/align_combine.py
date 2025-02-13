#!/usr/bin/python3

import os, sys

def usage(cmdline): #pylint: disable=unused-variable
    
    cmdline.set_author('Diana L. Giraldo (diana.giraldofranco@uantwerpen.be)')
    cmdline.set_synopsis('Perform iterative alignment, upsampling and averaging of multiple images.')
    #cmdline.add_description('This script is intended to be a baseline interpolation approach for multi-image super-resolution')
    cmdline.add_description('It uses mrtransform, mrregister, for_each and maskfilter from mrtrix3, and hd-bet.')
    
    cmdline.add_argument('inputs', help='Input images', nargs='+')
    cmdline.add_argument('grid', help='Image with template grid')
    cmdline.add_argument('output', help='output image')
    cmdline.add_argument('-masks', help='Masks of input images', nargs='+')
    cmdline.add_argument('-interp', type=str, choices=["nearest", "linear", "cubic", "sinc"], default='cubic', help='Interpolation method. Default: cubic.')
    cmdline.add_argument('-iter', type=int, help='Number of iterations. Default: 3.', default=3)
    
def execute(): #pylint: disable=unused-variable
    from mrtrix3 import app, image, path, run #pylint: disable=no-name-in-module, import-outside-toplevel
    
    app.check_output_path(app.ARGS.output)
    
    nLR = len(app.ARGS.inputs)
    app.console('Number of input images: ' + str(nLR))
    
    # Make scratch directory 
    app.make_scratch_dir()
    app.goto_scratch_dir()
    
    path.make_dir('inputs')
    path.make_dir('fovs')
    path.make_dir('regrid_inputs')
    path.make_dir('regrid_fovs')
    
    # Copy data to scratch directory
    run.command('mrconvert ' + path.from_user(app.ARGS.grid) + ' grid.nii -config RealignTransform 0')
    # PENDING: print size and spacing
    
    for i,imgpath in enumerate(app.ARGS.inputs):
        run.command('mrconvert ' + path.from_user(imgpath) + ' ' + os.path.join('inputs','img' + str(i) + '.nii'))
        run.command('mrcalc ' + path.from_user(imgpath) + ' -isnan -not ' + os.path.join('fovs','img' + str(i) + '.nii'))

    use_masks = False
    if app.ARGS.masks is not None:
        use_masks = True
        path.make_dir('masks')
        for i,imgpath in enumerate(app.ARGS.masks):
            run.command('mrconvert ' + path.from_user(imgpath) + ' ' + os.path.join('masks','img' + str(i) + '.nii'))
        
    # Start iterations
    for it in range(0, app.ARGS.iter):
        app.console('Starting iteration ' + str(it+1))
        
        if (it < 1):
            ref_img = os.path.join('inputs','img0.nii')
            ref_mask = os.path.join('masks','img0.nii')
            grid_temp = 'grid.nii'

        else:
            ref_img = 'output.nii.gz' 
            ref_mask = 'bet_mask.nii.gz'
            grid_temp = ref_img
        
        path.make_dir('transforms' + str(it))
        
        for img in range(nLR):
            transform = os.path.join('transforms' + str(it), 'img' + str(img) + '.txt')
            opt_mov_mask = ' -mask1 ' + os.path.join('masks', 'img' + str(img) + '.nii') if use_masks else '' 
            opt_ref_mask = ' -mask2 ' + ref_mask if use_masks else ''
            opt_rig_init = '' if (it < 1) else  ' -rigid_init_matrix ' + os.path.join('transforms' + str(it-1), 'img' + str(img) + '.txt')
            # Register
            run.command('mrregister inputs/img' + str(img) + '.nii ' + ref_img + ' -type rigid -rigid ' + transform + opt_mov_mask + opt_ref_mask + opt_rig_init)
            # Transform and regrid
            run.command('mrtransform inputs/img' + str(img) + '.nii regrid_inputs/img' + str(img) + '.nii  -linear ' + transform + ' -template ' + grid_temp + ' -interp ' + app.ARGS.interp + ' -force')
            run.command('mrtransform fovs/img' + str(img) + '.nii regrid_fovs/img' + str(img) + '.nii  -linear ' + transform + ' -template ' + grid_temp + ' -interp nearest -force')
            
             
        # Average (weighted by FOV)
        run.command('mrcat ' + ' '.join(['regrid_inputs/img' + str(img) + '.nii' for img in range(nLR)]) + ' - | mrmath - sum tmp_sum.nii -axis 3 -force')
        run.command('mrcat ' + ' '.join(['regrid_fovs/img' + str(img) + '.nii' for img in range(nLR)]) + ' - | mrmath - sum tmp_fovsum.nii -axis 3 -force')
        run.command('mrcalc tmp_sum.nii tmp_fovsum.nii -div output.nii -force')
        run.command('mrcalc output.nii -finite output.nii 0 -if -abs output.nii.gz -force')
        
        # Mask output
        if it < (app.ARGS.iter-1) and use_masks:
            os.system('hd-bet -i output.nii.gz -o bet.nii.gz -device cpu -mode fast -tta 0 > /dev/null')
            
    # Create output
    run.command('mrconvert output.nii.gz ' + path.from_user(app.ARGS.output),
                force=app.FORCE_OVERWRITE)
            

# Execute the script
import mrtrix3  #pylint: disable=wrong-import-position
mrtrix3.execute() #pylint: disable=no-member
            


    
    
    
