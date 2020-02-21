"""
Tests for resistivity model to geotiff.
"""
import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv

SIRSAM_RF = 'sirsam_Na_randomforest'


@pytest.fixture
def model_grid():
    grid = [-60, -30, -20, -15, -12, -10., -9., -8., -7., -6., -5., -4., -3., -2., -1., 0.,
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12, 15, 20, 30, 60]
    pad = 5
    cell_size = 1.
    return grid, pad, cell_size


@pytest.fixture
def grid_centres():
    centres = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5,
               -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11., 13.5,
               17.5, 25., 45.]
    stripped = [-9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5,
                4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
    stripped_keep_start = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5,
                           -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
                           9.5]
    return centres, stripped, stripped_keep_start


def test_get_centre_points(model_grid, grid_centres):
    grid, _, _ = model_grid
    ce = conv._get_centres(grid)

    # Should have one less centre point then there are nodes as the
    #  last node is the terminating boundary and has no centre.
    np.testing.assert_array_equal(ce, grid_centres[0])


def test_strip_padding(model_grid, grid_centres):
    _, pad, _ = model_grid
    centres, stripped, stripped_keep_start = grid_centres

    s = conv._strip_padding(centres, pad)
    # Stripped array should have padding removed from either side
    np.testing.assert_array_equal(s, stripped)

    s = conv._strip_padding(centres, pad, keep_start=True)
    # Stripped array while keeping start (for stripping Z padding)
    #  should only have padding removed from end
    np.testing.assert_array_equal(s, stripped_keep_start)

#class TestShiftmap:
#    IMAGE_OUTPUTS = [ 
#        SIRSAM_RF + '_shiftmap_most_likely.tif', 
#        SIRSAM_RF + '_shiftmap_most_likely_thumbnail.tif',
#        SIRSAM_RF + '_shiftmap_query_0.tif',
#        SIRSAM_RF + '_shiftmap_query_0_thumbnail.tif',
#        SIRSAM_RF + '_shiftmap_training_1.tif',
#        SIRSAM_RF + '_shiftmap_training_1_thumbnail.tif'
#    ]
#
#    OTHER_OUTPUTS = [
#        SIRSAM_RF + '_shiftmap_generated_points.txt',
#    ]
#
#    @staticmethod
#    @pytest.fixture(scope='class', autouse=True)
#    def run_sirsam_random_forest_shiftmap(request, num_procs, num_parts, sirsam_rf_conf, sirsam_rf_out):
#        """
#        Run the top level 'learn' command'. Removes generated output on
#        completion.
#        """
#        def finalize():
#            if os.path.exists(sirsam_rf_out):
#                shutil.rmtree(sirsam_rf_out)
#
#        request.addfinalizer(finalize)
#    
#        # If running with one processor, call uncoverml directly
#        if num_procs == 1:
#            try:
#                uncoverml.shiftmap([sirsam_rf_conf, '-p', num_parts])
#            # Catch SystemExit that gets raised by Click on competion
#            except SystemExit:
#                pass   
#        else:
#            try:
#                cmd = ['mpirun', '-n', str(num_procs),
#                       'uncoverml', 'shiftmap', sirsam_rf_conf, '-p', str(num_parts)]
#                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            except subprocess.CalledProcessError as e:
#                raise RuntimeError(f"'{cmd}' failed with error {e.returncode}: {e.output}")
#
#    @staticmethod
#    @pytest.fixture(params=IMAGE_OUTPUTS + OTHER_OUTPUTS)
#    def sirsam_rf_output(request, sirsam_rf_out):
#        return os.path.join(sirsam_rf_out, request.param)
#
#
#    @staticmethod
#    def test_output_exists(sirsam_rf_output):
#        """
#        Test that excepted outputs of 'shiftmap' command exist after
#        running.
#        """
#        assert os.path.exists(sirsam_rf_output)
#
#    @staticmethod
#    @pytest.fixture(params=IMAGE_OUTPUTS)
#    def sirsam_rf_comp_outputs(request, sirsam_rf_out, sirsam_rf_precomp_shiftmap):
#        """
#        """
#        return (
#            os.path.join(sirsam_rf_out, request.param),
#            os.path.join(sirsam_rf_precomp_shiftmap, request.param)
#        )
#
#    @staticmethod
#    def test_outputs_equal(sirsam_rf_comp_outputs): 
#        test = sirsam_rf_comp_outputs[0]
#        ref = sirsam_rf_comp_outputs[1]
#        
#        with rasterio.open(test) as test_img, rasterio.open(ref) as ref_img:
#            test_img_ar = test_img.read()
#            ref_img_ar = ref_img.read()
#          
#        assert np.array_equal(test_img_ar, ref_img_ar)
#
#
#class TestLearnCommand:
#    """
#    Tests the 'learn' command of the UncomverML CLI.
#    """
#    # Group the outputs of the learn command by filetype to make them easier to test.
#    SIRSAM_RF_JSON_OUTPUT = [
#        SIRSAM_RF + '_crossval_scores.json',
#        SIRSAM_RF + '_featureranks.json'
#    ]
#
#    SIRSAM_RF_FEATURE_DATA = 'features.pk'
#    SIRSAM_RF_TARGET_DATA = 'targets.pk'
#    SIRSAM_RF_MODEL = SIRSAM_RF + '.model'
#
#    SIRSAM_RF_CSV_OUTPUT = [
#        SIRSAM_RF + '_crossval_results.csv',
#        SIRSAM_RF + '_rawcovariates.csv',
#        SIRSAM_RF + '_rawcovariates_mask.csv'
#    ]
#
#    SIRSAM_RF_IMAGE_OUTPUT = [
#        SIRSAM_RF + '_featureranks.png',
#        SIRSAM_RF + '_featurerank_curves.png',
#        SIRSAM_RF + '_intersected.png',
#        SIRSAM_RF + '_correlation.png',
#        SIRSAM_RF + '_target_scaling.png',
#        SIRSAM_RF + '_real_vs_pred.png',
#        SIRSAM_RF + '_residual.png'
#    ]
#     
#   
#    SIRSAM_RF_OUTPUTS = [SIRSAM_RF_FEATURE_DATA, SIRSAM_RF_TARGET_DATA, 
#                         SIRSAM_RF_MODEL]
#    SIRSAM_RF_OUTPUTS += SIRSAM_RF_CSV_OUTPUT + SIRSAM_RF_IMAGE_OUTPUT + SIRSAM_RF_JSON_OUTPUT
#
#    @staticmethod
#    @pytest.fixture(scope='class', autouse=True)
#    def run_sirsam_random_forest_learning(request, num_procs, num_parts, sirsam_rf_conf, sirsam_rf_out):
#        """
#        Run the top level 'learn' command'. Removes generated output on
#        completion.
#        """
#        def finalize():
#            if os.path.exists(sirsam_rf_out):
#                shutil.rmtree(sirsam_rf_out)
#
#        request.addfinalizer(finalize)
#    
#        # If running with one processor, call uncoverml directly
#        if num_procs == 1:
#            try:
#                uncoverml.learn([sirsam_rf_conf, '-p', num_parts])
#            # Catch SystemExit that gets raised by Click on competion
#            except SystemExit:
#                pass   
#        else:
#            try:
#                cmd = ['mpirun', '-n', str(num_procs),
#                       'uncoverml', 'learn', sirsam_rf_conf, '-p', str(num_parts)]
#                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            except subprocess.CalledProcessError as e:
#                raise RuntimeError(f"'{cmd}' failed with error {e.returncode}: {e.output}")
#    
#    @staticmethod
#    @pytest.fixture(params=SIRSAM_RF_OUTPUTS)
#    def sirsam_rf_output(request, sirsam_rf_out):
#        return os.path.join(sirsam_rf_out, request.param)
#    
#    @staticmethod
#    def test_output_exists(sirsam_rf_output):
#        """
#        Test that excepted outputs of 'learn' command exist after
#        running.
#        """
#        assert os.path.exists(sirsam_rf_output)
#
#    @staticmethod
#    @pytest.fixture(params=SIRSAM_RF_CSV_OUTPUT)
#    def sirsam_rf_csv_outputs(request, sirsam_rf_out, sirsam_rf_precomp_learn):
#        return (
#            os.path.join(sirsam_rf_out, request.param),
#            os.path.join(sirsam_rf_precomp_learn,  request.param)
#        )
#
#    @staticmethod
#    def test_csv_outputs_match(sirsam_rf_csv_outputs):
#        """
#        Test that CSV covariate info matches precomputed output.
#        """
#        with open(sirsam_rf_csv_outputs[0]) as test, \
#                open(sirsam_rf_csv_outputs[1]) as precomp:
#            tl = test.readlines()
#            pl = precomp.readlines()
#        assert tl == pl
#
#    @staticmethod
#    @pytest.fixture(params=SIRSAM_RF_JSON_OUTPUT)
#    def sirsam_rf_json_outputs(request, sirsam_rf_out, sirsam_rf_precomp_learn):
#        return (
#            os.path.join(sirsam_rf_out, request.param),
#            os.path.join(sirsam_rf_precomp_learn, request.param)
#        )
#
#    @staticmethod
#    def test_json_outputs_match(sirsam_rf_json_outputs):
#        """
#        Test that JSON scores matches precomputed output.
#        """
#        with open(sirsam_rf_json_outputs[0]) as tf, open(sirsam_rf_json_outputs[1]) as pf:
#            test_json = json.load(tf)
#            precomp_json = json.load(pf)
#        assert test_json == precomp_json
#
#    @classmethod
#    def test_model_outputs_match(cls, sirsam_rf_out, sirsam_rf_precomp_learn):
#        """
#        Test that generated model matches precomputed model.
#        """
#        t_dict = _unpickle(os.path.join(sirsam_rf_out, cls.SIRSAM_RF_MODEL))
#        p_dict = _unpickle(os.path.join(sirsam_rf_precomp_learn, cls.SIRSAM_RF_MODEL))
#        assert t_dict['model'] == p_dict['model']
#        assert t_dict['config'] == p_dict['config']
#
#    @classmethod
#    def test_training_data_matches(cls, sirsam_rf_out, sirsam_rf_precomp_learn):
#        """
#        Test that pickled training data matches precomputed output.
#        """
#        t_image_chunk_sets, t_transform_sets, t_targets = \
#            _unpickle(os.path.join(sirsam_rf_out, cls.SIRSAM_RF_TRAINING_DATA))
#        p_image_chunk_sets, p_transform_sets, p_targets = \
#            _unpickle(os.path.join(sirsam_rf_precomp_learn, cls.SIRSAM_RF_TRAINING_DATA))
#        assert t_image_chunk_sets == p_image_chunk_sets
#        assert t_transform_sets == p_transform_sets
#        assert t_targets == p_targets
#
#    @classmethod
#    def test_pickled_features_match(cls, sirsam_rf_out, sirsam_rf_precomp_learn):
#        """
#        Test that pickled features match precomputed output.
#        """
#        t_features = _unpickle(os.path.join(sirsam_rf_out, cls.SIRSAM_RF_FEATURE_DATA))
#        p_features = _unpickle(os.path.join(sirsam_rf_precomp_learn, cls.SIRSAM_RF_FEATURE_DATA))
#        assert np.array_equal(t_features, p_features)
#
#    @classmethod
#    def test_pickled_targets_match(cls, sirsam_rf_out, sirsam_rf_precomp_learn):
#        """
#        Test that pickled targets match precomputed output.
#        """
#        t_targets = _unpickle(os.path.join(sirsam_rf_out, cls.SIRSAM_RF_TARGET_DATA))
#        p_targets = _unpickle(os.path.join(sirsam_rf_precomp_learn, cls.SIRSAM_RF_TARGET_DATA))
#        assert t_targets == p_targets
#
#    @classmethod
#    def test_multi_random_forest_caching(cls, sirsam_rf_out):
#        model = _unpickle(os.path.join(sirsam_rf_out, cls.SIRSAM_RF_MODEL))
#        assert model._randomforests == model.n_estimators
#
#def _unpickle(path):
#    with open(path, 'rb') as f:
#        return pickle.load(f)
#
#@pytest.mark.predict
#class TestPredictCommand:
#    SIRSAM_PREDICTION_MAPS = [
#        SIRSAM_RF + '_lower_quantile.tif',
#        SIRSAM_RF + '_lower_quantile_thumbnail.tif',
#        SIRSAM_RF + '_prediction.tif',
#        SIRSAM_RF + '_prediction_thumbnail.tif',
#        SIRSAM_RF + '_upper_quantile.tif',
#        SIRSAM_RF + '_upper_quantile_thumbnail.tif',
#        SIRSAM_RF + '_variance.tif',
#        SIRSAM_RF + '_variance_thumbnail.tif'
#    ]
#    
#    SIRSAM_RF_IMAGE_OUTPUT = [
#        #SIRSAM_RF + '_real_vs_pred.png',
#        #SIRSAM_RF + '_residual.png'
#    ]
#        
#    SIRSAM_RF_MF_METADATA = SIRSAM_RF + '_metadata.txt'
#
#    @staticmethod
#    @pytest.fixture(scope='class', autouse=True)
#    def run_sirsam_random_forest_prediction(request, num_procs, num_parts, sirsam_rf_out, sirsam_rf_conf, 
#                                            sirsam_rf_precomp_learn):
#        """
#        Run the top level 'predict' command'. Removes generated output on
#        completion.
#        """
#        def finalize():
#            if os.path.exists(sirsam_rf_out):
#                shutil.rmtree(sirsam_rf_out)
#
#        request.addfinalizer(finalize)
#
#        # Copy precomputed files from learn step to the output directory
#        shutil.copytree(sirsam_rf_precomp_learn, sirsam_rf_out)
#
#        # If running with one processor, call uncoverml directly
#        if num_procs == 1:
#            try:
#                uncoverml.predict([sirsam_rf_conf, '-p', num_parts])
#            # Catch SystemExit that gets raised by Click on competion
#            except SystemExit:
#                pass   
#        else:
#            try:
#                cmd = ['mpirun', '-n', str(num_procs),
#                       'uncoverml', 'predict', sirsam_rf_conf, '-p', str(num_parts)]
#                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#            except subprocess.CalledProcessError as e:
#                raise RuntimeError(f"'{cmd}' failed with error {e.returncode}: {e.output}")
#
#    @staticmethod
#    @pytest.fixture(params=SIRSAM_PREDICTION_MAPS + [SIRSAM_RF_MF_METADATA] + SIRSAM_RF_IMAGE_OUTPUT)
#    def sirsam_rf_output(request, sirsam_rf_out):
#        return os.path.join(sirsam_rf_out, request.param)
#    
#    @staticmethod
#    def test_output_exists(sirsam_rf_output):
#        """
#        Test that excepted outputs of 'predict' command exist after
#        running.
#        """
#        assert os.path.exists(sirsam_rf_output)

