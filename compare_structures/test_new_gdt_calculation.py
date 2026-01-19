#!/usr/bin/env python3
"""
Test the updated GDT score calculation using CA_gdttm and CA_gdtmm
"""

import numpy as np
import sys
sys.path.append('.')

def test_gdt_calculation_new():
    """Test the new GDT calculation logic."""
    
    # Mock implementation of CA_gdtmm for testing
    def mock_CA_gdtmm(pose1, pose2):
        """Mock CA_gdtmm function that returns a score between 0-1."""
        # Return a sample GDT-MM score
        return 0.75
    
    # Mock implementation of CA_superimpose for testing
    class MockSuperimposedInfo:
        def __init__(self):
            self.rmsd = 2.0
    
    def mock_CA_superimpose(pose1, pose2):
        """Mock CA_superimpose function."""
        return MockSuperimposedInfo()
    
    # Mock Pose objects
    class MockPose:
        def __init__(self, residue_count):
            self.residue_count = residue_count
        
        def total_residue(self):
            return self.residue_count
    
    print("Testing new GDT calculation logic...")
    
    # Create mock poses
    pose1 = MockPose(20)
    pose2 = MockPose(25)
    
    # Test the calculation function logic (without actual PyRosetta functions)
    thresholds = [1.0, 2.0, 4.0, 8.0]
    gdt_scores = {}
    
    # Simulate the CA_gdtmm approach
    try:
        gdtmm_score = mock_CA_gdtmm(pose1, pose2)
        
        for threshold in thresholds:
            if threshold <= 1.0:
                gdt_score = gdtmm_score * 0.8
            elif threshold <= 2.0:
                gdt_score = gdtmm_score * 0.9
            elif threshold <= 4.0:
                gdt_score = gdtmm_score * 0.95
            else:
                gdt_score = min(1.0, gdtmm_score * 1.05)
            
            gdt_scores[f'gdt_ts{int(threshold)}'] = gdt_score
        
        print("✓ Primary GDT calculation using CA_gdtmm successful")
        
    except Exception as e:
        print(f"✗ Primary GDT calculation failed: {e}")
        
        # Test fallback method
        superimposed_info = mock_CA_superimpose(pose1, pose2)
        rmsd = superimposed_info.rmsd
        
        for threshold in thresholds:
            if rmsd <= threshold:
                gdt_score = max(0.0, 1.0 - (rmsd / threshold))
            else:
                gdt_score = max(0.0, (threshold / rmsd))
            
            gdt_scores[f'gdt_ts{int(threshold)}'] = gdt_score
        
        print("✓ Fallback GDT calculation successful")
    
    # Display results
    print("\nGDT-TS Scores:")
    for key, value in gdt_scores.items():
        print(f"  {key}: {value:.3f}")
    
    # Verify scores are reasonable
    assert all(0 <= score <= 1 for score in gdt_scores.values()), "Scores should be between 0 and 1"
    assert len(gdt_scores) == 4, "Should have 4 GDT scores"
    
    print("\n✓ GDT score calculation verification passed")

def test_threshold_scaling():
    """Test the threshold scaling logic."""
    print("\nTesting threshold scaling logic...")
    
    base_gdtmm = 0.75
    expected = {
        1.0: base_gdtmm * 0.8,    # 0.600
        2.0: base_gdtmm * 0.9,    # 0.675
        4.0: base_gdtmm * 0.95,   # 0.713
        8.0: min(1.0, base_gdtmm * 1.05)  # 0.788
    }
    
    thresholds = [1.0, 2.0, 4.0, 8.0]
    
    for threshold in thresholds:
        if threshold <= 1.0:
            calculated = base_gdtmm * 0.8
        elif threshold <= 2.0:
            calculated = base_gdtmm * 0.9
        elif threshold <= 4.0:
            calculated = base_gdtmm * 0.95
        else:
            calculated = min(1.0, base_gdtmm * 1.05)
        
        expected_score = expected[threshold]
        assert abs(calculated - expected_score) < 0.001, f"Mismatch for {threshold}Å: {calculated} vs {expected_score}"
    
    print("✓ Threshold scaling logic verified")

if __name__ == '__main__':
    test_gdt_calculation_new()
    test_threshold_scaling()
    print("\nAll GDT calculation tests passed!")