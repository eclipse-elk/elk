/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData.impl;

import org.eclipse.elk.core.meta.metaData.MdAlgorithm;
import org.eclipse.elk.core.meta.metaData.MdBundle;
import org.eclipse.elk.core.meta.metaData.MdBundleMember;
import org.eclipse.elk.core.meta.metaData.MdCategory;
import org.eclipse.elk.core.meta.metaData.MdGraphFeature;
import org.eclipse.elk.core.meta.metaData.MdGroup;
import org.eclipse.elk.core.meta.metaData.MdGroupOrOption;
import org.eclipse.elk.core.meta.metaData.MdModel;
import org.eclipse.elk.core.meta.metaData.MdOption;
import org.eclipse.elk.core.meta.metaData.MdOptionDependency;
import org.eclipse.elk.core.meta.metaData.MdOptionSupport;
import org.eclipse.elk.core.meta.metaData.MdOptionTargetType;
import org.eclipse.elk.core.meta.metaData.MetaDataFactory;
import org.eclipse.elk.core.meta.metaData.MetaDataPackage;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EEnum;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.emf.ecore.EReference;

import org.eclipse.emf.ecore.impl.EPackageImpl;

import org.eclipse.xtext.common.types.TypesPackage;

import org.eclipse.xtext.xbase.XbasePackage;

import org.eclipse.xtext.xtype.XtypePackage;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model <b>Package</b>.
 * <!-- end-user-doc -->
 * @generated
 */
public class MetaDataPackageImpl extends EPackageImpl implements MetaDataPackage
{
  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdModelEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdBundleEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdBundleMemberEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdGroupOrOptionEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdGroupEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdOptionEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdOptionDependencyEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdAlgorithmEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdCategoryEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EClass mdOptionSupportEClass = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum mdOptionTargetTypeEEnum = null;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private EEnum mdGraphFeatureEEnum = null;

  /**
   * Creates an instance of the model <b>Package</b>, registered with
   * {@link org.eclipse.emf.ecore.EPackage.Registry EPackage.Registry} by the package
   * package URI value.
   * <p>Note: the correct way to create the package is via the static
   * factory method {@link #init init()}, which also performs
   * initialization of the package, or returns the registered package,
   * if one already exists.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.emf.ecore.EPackage.Registry
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#eNS_URI
   * @see #init()
   * @generated
   */
  private MetaDataPackageImpl()
  {
    super(eNS_URI, MetaDataFactory.eINSTANCE);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private static boolean isInited = false;

  /**
   * Creates, registers, and initializes the <b>Package</b> for this model, and for any others upon which it depends.
   *
   * <p>This method is used to initialize {@link MetaDataPackage#eINSTANCE} when that field is accessed.
   * Clients should not invoke it directly. Instead, they should simply access that field to obtain the package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #eNS_URI
   * @see #createPackageContents()
   * @see #initializePackageContents()
   * @generated
   */
  public static MetaDataPackage init()
  {
    if (isInited) return (MetaDataPackage)EPackage.Registry.INSTANCE.getEPackage(MetaDataPackage.eNS_URI);

    // Obtain or create and register package
    Object registeredMetaDataPackage = EPackage.Registry.INSTANCE.get(eNS_URI);
    MetaDataPackageImpl theMetaDataPackage = registeredMetaDataPackage instanceof MetaDataPackageImpl ? (MetaDataPackageImpl)registeredMetaDataPackage : new MetaDataPackageImpl();

    isInited = true;

    // Initialize simple dependencies
    TypesPackage.eINSTANCE.eClass();
    XbasePackage.eINSTANCE.eClass();
    XtypePackage.eINSTANCE.eClass();

    // Create package meta-data objects
    theMetaDataPackage.createPackageContents();

    // Initialize created meta-data
    theMetaDataPackage.initializePackageContents();

    // Mark meta-data to indicate it can't be changed
    theMetaDataPackage.freeze();

    // Update the registry and return the package
    EPackage.Registry.INSTANCE.put(MetaDataPackage.eNS_URI, theMetaDataPackage);
    return theMetaDataPackage;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdModel()
  {
    return mdModelEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdModel_Name()
  {
    return (EAttribute)mdModelEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdModel_ImportSection()
  {
    return (EReference)mdModelEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdModel_Bundle()
  {
    return (EReference)mdModelEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdBundle()
  {
    return mdBundleEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundle_Label()
  {
    return (EAttribute)mdBundleEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundle_TargetClass()
  {
    return (EAttribute)mdBundleEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundle_DocumentationFolder()
  {
    return (EAttribute)mdBundleEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundle_IdPrefix()
  {
    return (EAttribute)mdBundleEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdBundle_Members()
  {
    return (EReference)mdBundleEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdBundleMember()
  {
    return mdBundleMemberEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundleMember_Name()
  {
    return (EAttribute)mdBundleMemberEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdBundleMember_Documentation()
  {
    return (EAttribute)mdBundleMemberEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdGroupOrOption()
  {
    return mdGroupOrOptionEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdGroup()
  {
    return mdGroupEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdGroup_Children()
  {
    return (EReference)mdGroupEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdOption()
  {
    return mdOptionEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Deprecated()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Advanced()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Programmatic()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Output()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Global()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOption_Type()
  {
    return (EReference)mdOptionEClass.getEStructuralFeatures().get(5);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Label()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(6);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Description()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(7);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOption_DefaultValue()
  {
    return (EReference)mdOptionEClass.getEStructuralFeatures().get(8);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOption_LowerBound()
  {
    return (EReference)mdOptionEClass.getEStructuralFeatures().get(9);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOption_UpperBound()
  {
    return (EReference)mdOptionEClass.getEStructuralFeatures().get(10);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_Targets()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(11);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOption_LegacyIds()
  {
    return (EAttribute)mdOptionEClass.getEStructuralFeatures().get(12);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOption_Dependencies()
  {
    return (EReference)mdOptionEClass.getEStructuralFeatures().get(13);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdOptionDependency()
  {
    return mdOptionDependencyEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOptionDependency_Target()
  {
    return (EReference)mdOptionDependencyEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOptionDependency_Value()
  {
    return (EReference)mdOptionDependencyEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdAlgorithm()
  {
    return mdAlgorithmEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_Deprecated()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdAlgorithm_Provider()
  {
    return (EReference)mdAlgorithmEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_Parameter()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_Label()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(3);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_TargetClass()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(4);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_Description()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(5);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdAlgorithm_Category()
  {
    return (EReference)mdAlgorithmEClass.getEStructuralFeatures().get(6);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_PreviewImage()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(7);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdAlgorithm_SupportedFeatures()
  {
    return (EAttribute)mdAlgorithmEClass.getEStructuralFeatures().get(8);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdAlgorithm_Validator()
  {
    return (EReference)mdAlgorithmEClass.getEStructuralFeatures().get(9);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdAlgorithm_SupportedOptions()
  {
    return (EReference)mdAlgorithmEClass.getEStructuralFeatures().get(10);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdCategory()
  {
    return mdCategoryEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdCategory_Deprecated()
  {
    return (EAttribute)mdCategoryEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdCategory_Label()
  {
    return (EAttribute)mdCategoryEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdCategory_Description()
  {
    return (EAttribute)mdCategoryEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EClass getMdOptionSupport()
  {
    return mdOptionSupportEClass;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOptionSupport_Option()
  {
    return (EReference)mdOptionSupportEClass.getEStructuralFeatures().get(0);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EReference getMdOptionSupport_Value()
  {
    return (EReference)mdOptionSupportEClass.getEStructuralFeatures().get(1);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EAttribute getMdOptionSupport_Documentation()
  {
    return (EAttribute)mdOptionSupportEClass.getEStructuralFeatures().get(2);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getMdOptionTargetType()
  {
    return mdOptionTargetTypeEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EEnum getMdGraphFeature()
  {
    return mdGraphFeatureEEnum;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MetaDataFactory getMetaDataFactory()
  {
    return (MetaDataFactory)getEFactoryInstance();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private boolean isCreated = false;

  /**
   * Creates the meta-model objects for the package.  This method is
   * guarded to have no affect on any invocation but its first.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void createPackageContents()
  {
    if (isCreated) return;
    isCreated = true;

    // Create classes and their features
    mdModelEClass = createEClass(MD_MODEL);
    createEAttribute(mdModelEClass, MD_MODEL__NAME);
    createEReference(mdModelEClass, MD_MODEL__IMPORT_SECTION);
    createEReference(mdModelEClass, MD_MODEL__BUNDLE);

    mdBundleEClass = createEClass(MD_BUNDLE);
    createEAttribute(mdBundleEClass, MD_BUNDLE__LABEL);
    createEAttribute(mdBundleEClass, MD_BUNDLE__TARGET_CLASS);
    createEAttribute(mdBundleEClass, MD_BUNDLE__DOCUMENTATION_FOLDER);
    createEAttribute(mdBundleEClass, MD_BUNDLE__ID_PREFIX);
    createEReference(mdBundleEClass, MD_BUNDLE__MEMBERS);

    mdBundleMemberEClass = createEClass(MD_BUNDLE_MEMBER);
    createEAttribute(mdBundleMemberEClass, MD_BUNDLE_MEMBER__NAME);
    createEAttribute(mdBundleMemberEClass, MD_BUNDLE_MEMBER__DOCUMENTATION);

    mdGroupOrOptionEClass = createEClass(MD_GROUP_OR_OPTION);

    mdGroupEClass = createEClass(MD_GROUP);
    createEReference(mdGroupEClass, MD_GROUP__CHILDREN);

    mdOptionEClass = createEClass(MD_OPTION);
    createEAttribute(mdOptionEClass, MD_OPTION__DEPRECATED);
    createEAttribute(mdOptionEClass, MD_OPTION__ADVANCED);
    createEAttribute(mdOptionEClass, MD_OPTION__PROGRAMMATIC);
    createEAttribute(mdOptionEClass, MD_OPTION__OUTPUT);
    createEAttribute(mdOptionEClass, MD_OPTION__GLOBAL);
    createEReference(mdOptionEClass, MD_OPTION__TYPE);
    createEAttribute(mdOptionEClass, MD_OPTION__LABEL);
    createEAttribute(mdOptionEClass, MD_OPTION__DESCRIPTION);
    createEReference(mdOptionEClass, MD_OPTION__DEFAULT_VALUE);
    createEReference(mdOptionEClass, MD_OPTION__LOWER_BOUND);
    createEReference(mdOptionEClass, MD_OPTION__UPPER_BOUND);
    createEAttribute(mdOptionEClass, MD_OPTION__TARGETS);
    createEAttribute(mdOptionEClass, MD_OPTION__LEGACY_IDS);
    createEReference(mdOptionEClass, MD_OPTION__DEPENDENCIES);

    mdOptionDependencyEClass = createEClass(MD_OPTION_DEPENDENCY);
    createEReference(mdOptionDependencyEClass, MD_OPTION_DEPENDENCY__TARGET);
    createEReference(mdOptionDependencyEClass, MD_OPTION_DEPENDENCY__VALUE);

    mdAlgorithmEClass = createEClass(MD_ALGORITHM);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__DEPRECATED);
    createEReference(mdAlgorithmEClass, MD_ALGORITHM__PROVIDER);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__PARAMETER);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__LABEL);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__TARGET_CLASS);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__DESCRIPTION);
    createEReference(mdAlgorithmEClass, MD_ALGORITHM__CATEGORY);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__PREVIEW_IMAGE);
    createEAttribute(mdAlgorithmEClass, MD_ALGORITHM__SUPPORTED_FEATURES);
    createEReference(mdAlgorithmEClass, MD_ALGORITHM__VALIDATOR);
    createEReference(mdAlgorithmEClass, MD_ALGORITHM__SUPPORTED_OPTIONS);

    mdCategoryEClass = createEClass(MD_CATEGORY);
    createEAttribute(mdCategoryEClass, MD_CATEGORY__DEPRECATED);
    createEAttribute(mdCategoryEClass, MD_CATEGORY__LABEL);
    createEAttribute(mdCategoryEClass, MD_CATEGORY__DESCRIPTION);

    mdOptionSupportEClass = createEClass(MD_OPTION_SUPPORT);
    createEReference(mdOptionSupportEClass, MD_OPTION_SUPPORT__OPTION);
    createEReference(mdOptionSupportEClass, MD_OPTION_SUPPORT__VALUE);
    createEAttribute(mdOptionSupportEClass, MD_OPTION_SUPPORT__DOCUMENTATION);

    // Create enums
    mdOptionTargetTypeEEnum = createEEnum(MD_OPTION_TARGET_TYPE);
    mdGraphFeatureEEnum = createEEnum(MD_GRAPH_FEATURE);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private boolean isInitialized = false;

  /**
   * Complete the initialization of the package and its meta-model.  This
   * method is guarded to have no affect on any invocation but its first.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void initializePackageContents()
  {
    if (isInitialized) return;
    isInitialized = true;

    // Initialize package
    setName(eNAME);
    setNsPrefix(eNS_PREFIX);
    setNsURI(eNS_URI);

    // Obtain other dependent packages
    XtypePackage theXtypePackage = (XtypePackage)EPackage.Registry.INSTANCE.getEPackage(XtypePackage.eNS_URI);
    TypesPackage theTypesPackage = (TypesPackage)EPackage.Registry.INSTANCE.getEPackage(TypesPackage.eNS_URI);
    XbasePackage theXbasePackage = (XbasePackage)EPackage.Registry.INSTANCE.getEPackage(XbasePackage.eNS_URI);

    // Create type parameters

    // Set bounds for type parameters

    // Add supertypes to classes
    mdGroupOrOptionEClass.getESuperTypes().add(this.getMdBundleMember());
    mdGroupEClass.getESuperTypes().add(this.getMdGroupOrOption());
    mdOptionEClass.getESuperTypes().add(this.getMdGroupOrOption());
    mdAlgorithmEClass.getESuperTypes().add(this.getMdBundleMember());
    mdCategoryEClass.getESuperTypes().add(this.getMdBundleMember());

    // Initialize classes and features; add operations and parameters
    initEClass(mdModelEClass, MdModel.class, "MdModel", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdModel_Name(), ecorePackage.getEString(), "name", null, 0, 1, MdModel.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdModel_ImportSection(), theXtypePackage.getXImportSection(), null, "importSection", null, 0, 1, MdModel.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdModel_Bundle(), this.getMdBundle(), null, "bundle", null, 0, 1, MdModel.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdBundleEClass, MdBundle.class, "MdBundle", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdBundle_Label(), ecorePackage.getEString(), "label", null, 0, 1, MdBundle.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdBundle_TargetClass(), ecorePackage.getEString(), "targetClass", null, 0, 1, MdBundle.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdBundle_DocumentationFolder(), ecorePackage.getEString(), "documentationFolder", null, 0, 1, MdBundle.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdBundle_IdPrefix(), ecorePackage.getEString(), "idPrefix", null, 0, 1, MdBundle.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdBundle_Members(), this.getMdBundleMember(), null, "members", null, 0, -1, MdBundle.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdBundleMemberEClass, MdBundleMember.class, "MdBundleMember", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdBundleMember_Name(), ecorePackage.getEString(), "name", null, 0, 1, MdBundleMember.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdBundleMember_Documentation(), ecorePackage.getEString(), "documentation", null, 0, 1, MdBundleMember.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdGroupOrOptionEClass, MdGroupOrOption.class, "MdGroupOrOption", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);

    initEClass(mdGroupEClass, MdGroup.class, "MdGroup", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getMdGroup_Children(), this.getMdGroupOrOption(), null, "children", null, 0, -1, MdGroup.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdOptionEClass, MdOption.class, "MdOption", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdOption_Deprecated(), ecorePackage.getEBoolean(), "deprecated", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Advanced(), ecorePackage.getEBoolean(), "advanced", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Programmatic(), ecorePackage.getEBoolean(), "programmatic", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Output(), ecorePackage.getEBoolean(), "output", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Global(), ecorePackage.getEBoolean(), "global", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOption_Type(), theTypesPackage.getJvmTypeReference(), null, "type", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Label(), ecorePackage.getEString(), "label", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Description(), ecorePackage.getEString(), "description", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOption_DefaultValue(), theXbasePackage.getXExpression(), null, "defaultValue", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOption_LowerBound(), theXbasePackage.getXExpression(), null, "lowerBound", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOption_UpperBound(), theXbasePackage.getXExpression(), null, "upperBound", null, 0, 1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_Targets(), this.getMdOptionTargetType(), "targets", null, 0, -1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, !IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOption_LegacyIds(), ecorePackage.getEString(), "legacyIds", null, 0, -1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, !IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOption_Dependencies(), this.getMdOptionDependency(), null, "dependencies", null, 0, -1, MdOption.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdOptionDependencyEClass, MdOptionDependency.class, "MdOptionDependency", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getMdOptionDependency_Target(), this.getMdOption(), null, "target", null, 0, 1, MdOptionDependency.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOptionDependency_Value(), theXbasePackage.getXExpression(), null, "value", null, 0, 1, MdOptionDependency.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdAlgorithmEClass, MdAlgorithm.class, "MdAlgorithm", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdAlgorithm_Deprecated(), ecorePackage.getEBoolean(), "deprecated", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdAlgorithm_Provider(), theTypesPackage.getJvmTypeReference(), null, "provider", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_Parameter(), ecorePackage.getEString(), "parameter", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_Label(), ecorePackage.getEString(), "label", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_TargetClass(), ecorePackage.getEString(), "targetClass", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_Description(), ecorePackage.getEString(), "description", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdAlgorithm_Category(), this.getMdCategory(), null, "category", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_PreviewImage(), ecorePackage.getEString(), "previewImage", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdAlgorithm_SupportedFeatures(), this.getMdGraphFeature(), "supportedFeatures", null, 0, -1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, !IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdAlgorithm_Validator(), theTypesPackage.getJvmTypeReference(), null, "validator", null, 0, 1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdAlgorithm_SupportedOptions(), this.getMdOptionSupport(), null, "supportedOptions", null, 0, -1, MdAlgorithm.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdCategoryEClass, MdCategory.class, "MdCategory", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEAttribute(getMdCategory_Deprecated(), ecorePackage.getEBoolean(), "deprecated", null, 0, 1, MdCategory.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdCategory_Label(), ecorePackage.getEString(), "label", null, 0, 1, MdCategory.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdCategory_Description(), ecorePackage.getEString(), "description", null, 0, 1, MdCategory.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    initEClass(mdOptionSupportEClass, MdOptionSupport.class, "MdOptionSupport", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
    initEReference(getMdOptionSupport_Option(), this.getMdOption(), null, "option", null, 0, 1, MdOptionSupport.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEReference(getMdOptionSupport_Value(), theXbasePackage.getXExpression(), null, "value", null, 0, 1, MdOptionSupport.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
    initEAttribute(getMdOptionSupport_Documentation(), ecorePackage.getEString(), "documentation", null, 0, 1, MdOptionSupport.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

    // Initialize enums and add enum literals
    initEEnum(mdOptionTargetTypeEEnum, MdOptionTargetType.class, "MdOptionTargetType");
    addEEnumLiteral(mdOptionTargetTypeEEnum, MdOptionTargetType.PARENTS);
    addEEnumLiteral(mdOptionTargetTypeEEnum, MdOptionTargetType.NODES);
    addEEnumLiteral(mdOptionTargetTypeEEnum, MdOptionTargetType.EDGES);
    addEEnumLiteral(mdOptionTargetTypeEEnum, MdOptionTargetType.PORTS);
    addEEnumLiteral(mdOptionTargetTypeEEnum, MdOptionTargetType.LABELS);

    initEEnum(mdGraphFeatureEEnum, MdGraphFeature.class, "MdGraphFeature");
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.SELF_LOOPS);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.INSIDE_SELF_LOOPS);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.MULTI_EDGES);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.EDGE_LABELS);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.PORTS);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.COMPOUND);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.CLUSTERS);
    addEEnumLiteral(mdGraphFeatureEEnum, MdGraphFeature.DISCONNECTED);

    // Create resource
    createResource(eNS_URI);
  }

} //MetaDataPackageImpl
