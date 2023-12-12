from setuptools import setup, find_packages

setup(
    name='BS_predictor',
    version='2.0.0',
    description = "Machine Learning-Based Binding Site Prediction for Protein Structures",
    
    author="Denis Expósito Navarro, Marta Alonso Caubilla, Yolanda Andrés López",
    author_email="deexna4@gmail.com, marta.alonso07@estudiant.upf.edu, yolanda.andres01@estudiant.upf.edu",
    
    packages=find_packages(),
    
    package_data = {
        "BS_predictor": ["data/**/*"],
    },
    include_package_data=True,
    
    entry_points={
        'console_scripts': [
            'BS_predictor=BS_predictor.BS_predictor:main',
        ]
    }
    )

