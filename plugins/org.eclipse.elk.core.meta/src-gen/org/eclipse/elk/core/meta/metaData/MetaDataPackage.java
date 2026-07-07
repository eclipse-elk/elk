/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EEnum;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.emf.ecore.EReference;

/**
 * <!-- begin-user-doc -->
 * The <b>Package</b> for the model.
 * It contains accessors for the meta objects to represent
 * <ul>
 *   <li>each class,</li>
 *   <li>each feature of each class,</li>
 *   <li>each enum,</li>
 *   <li>and each data type</li>
 * </ul>
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.meta.metaData.MetaDataFactory
 * @model kind="package"
 * @generated
 */
public interface MetaDataPackage extends EPackage
{
  /**
   * The package name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNAME = "metaData";

  /**
   * The package namespace URI.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNS_URI = "http://www.eclipse.org/elk/MetaData";

  /**
   * The package namespace name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  String eNS_PREFIX = "metaData";

  /**
   * The singleton instance of the package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  MetaDataPackage eINSTANCE = org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl.init();

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdModelImpl <em>Md Model</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdModelImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdModel()
   * @generated
   */
  int MD_MODEL = 0;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_MODEL__NAME = 0;

  /**
   * The feature id for the '<em><b>Import Section</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_MODEL__IMPORT_SECTION = 1;

  /**
   * The feature id for the '<em><b>Bundle</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_MODEL__BUNDLE = 2;

  /**
   * The number of structural features of the '<em>Md Model</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_MODEL_FEATURE_COUNT = 3;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl <em>Md Bundle</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdBundle()
   * @generated
   */
  int MD_BUNDLE = 1;

  /**
   * The feature id for the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE__LABEL = 0;

  /**
   * The feature id for the '<em><b>Target Class</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE__TARGET_CLASS = 1;

  /**
   * The feature id for the '<em><b>Documentation Folder</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE__DOCUMENTATION_FOLDER = 2;

  /**
   * The feature id for the '<em><b>Id Prefix</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE__ID_PREFIX = 3;

  /**
   * The feature id for the '<em><b>Members</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE__MEMBERS = 4;

  /**
   * The number of structural features of the '<em>Md Bundle</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE_FEATURE_COUNT = 5;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleMemberImpl <em>Md Bundle Member</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdBundleMemberImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdBundleMember()
   * @generated
   */
  int MD_BUNDLE_MEMBER = 2;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE_MEMBER__NAME = 0;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE_MEMBER__DOCUMENTATION = 1;

  /**
   * The number of structural features of the '<em>Md Bundle Member</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_BUNDLE_MEMBER_FEATURE_COUNT = 2;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdGroupOrOptionImpl <em>Md Group Or Option</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdGroupOrOptionImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGroupOrOption()
   * @generated
   */
  int MD_GROUP_OR_OPTION = 3;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP_OR_OPTION__NAME = MD_BUNDLE_MEMBER__NAME;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP_OR_OPTION__DOCUMENTATION = MD_BUNDLE_MEMBER__DOCUMENTATION;

  /**
   * The number of structural features of the '<em>Md Group Or Option</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP_OR_OPTION_FEATURE_COUNT = MD_BUNDLE_MEMBER_FEATURE_COUNT + 0;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdGroupImpl <em>Md Group</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdGroupImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGroup()
   * @generated
   */
  int MD_GROUP = 4;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP__NAME = MD_GROUP_OR_OPTION__NAME;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP__DOCUMENTATION = MD_GROUP_OR_OPTION__DOCUMENTATION;

  /**
   * The feature id for the '<em><b>Children</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP__CHILDREN = MD_GROUP_OR_OPTION_FEATURE_COUNT + 0;

  /**
   * The number of structural features of the '<em>Md Group</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_GROUP_FEATURE_COUNT = MD_GROUP_OR_OPTION_FEATURE_COUNT + 1;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl <em>Md Option</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOption()
   * @generated
   */
  int MD_OPTION = 5;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__NAME = MD_GROUP_OR_OPTION__NAME;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__DOCUMENTATION = MD_GROUP_OR_OPTION__DOCUMENTATION;

  /**
   * The feature id for the '<em><b>Deprecated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__DEPRECATED = MD_GROUP_OR_OPTION_FEATURE_COUNT + 0;

  /**
   * The feature id for the '<em><b>Advanced</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__ADVANCED = MD_GROUP_OR_OPTION_FEATURE_COUNT + 1;

  /**
   * The feature id for the '<em><b>Programmatic</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__PROGRAMMATIC = MD_GROUP_OR_OPTION_FEATURE_COUNT + 2;

  /**
   * The feature id for the '<em><b>Output</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__OUTPUT = MD_GROUP_OR_OPTION_FEATURE_COUNT + 3;

  /**
   * The feature id for the '<em><b>Global</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__GLOBAL = MD_GROUP_OR_OPTION_FEATURE_COUNT + 4;

  /**
   * The feature id for the '<em><b>Type</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__TYPE = MD_GROUP_OR_OPTION_FEATURE_COUNT + 5;

  /**
   * The feature id for the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__LABEL = MD_GROUP_OR_OPTION_FEATURE_COUNT + 6;

  /**
   * The feature id for the '<em><b>Description</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__DESCRIPTION = MD_GROUP_OR_OPTION_FEATURE_COUNT + 7;

  /**
   * The feature id for the '<em><b>Default Value</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__DEFAULT_VALUE = MD_GROUP_OR_OPTION_FEATURE_COUNT + 8;

  /**
   * The feature id for the '<em><b>Lower Bound</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__LOWER_BOUND = MD_GROUP_OR_OPTION_FEATURE_COUNT + 9;

  /**
   * The feature id for the '<em><b>Upper Bound</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__UPPER_BOUND = MD_GROUP_OR_OPTION_FEATURE_COUNT + 10;

  /**
   * The feature id for the '<em><b>Targets</b></em>' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__TARGETS = MD_GROUP_OR_OPTION_FEATURE_COUNT + 11;

  /**
   * The feature id for the '<em><b>Legacy Ids</b></em>' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__LEGACY_IDS = MD_GROUP_OR_OPTION_FEATURE_COUNT + 12;

  /**
   * The feature id for the '<em><b>Dependencies</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION__DEPENDENCIES = MD_GROUP_OR_OPTION_FEATURE_COUNT + 13;

  /**
   * The number of structural features of the '<em>Md Option</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_FEATURE_COUNT = MD_GROUP_OR_OPTION_FEATURE_COUNT + 14;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionDependencyImpl <em>Md Option Dependency</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionDependencyImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionDependency()
   * @generated
   */
  int MD_OPTION_DEPENDENCY = 6;

  /**
   * The feature id for the '<em><b>Target</b></em>' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_DEPENDENCY__TARGET = 0;

  /**
   * The feature id for the '<em><b>Value</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_DEPENDENCY__VALUE = 1;

  /**
   * The number of structural features of the '<em>Md Option Dependency</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_DEPENDENCY_FEATURE_COUNT = 2;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl <em>Md Algorithm</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdAlgorithm()
   * @generated
   */
  int MD_ALGORITHM = 7;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__NAME = MD_BUNDLE_MEMBER__NAME;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__DOCUMENTATION = MD_BUNDLE_MEMBER__DOCUMENTATION;

  /**
   * The feature id for the '<em><b>Deprecated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__DEPRECATED = MD_BUNDLE_MEMBER_FEATURE_COUNT + 0;

  /**
   * The feature id for the '<em><b>Provider</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__PROVIDER = MD_BUNDLE_MEMBER_FEATURE_COUNT + 1;

  /**
   * The feature id for the '<em><b>Parameter</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__PARAMETER = MD_BUNDLE_MEMBER_FEATURE_COUNT + 2;

  /**
   * The feature id for the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__LABEL = MD_BUNDLE_MEMBER_FEATURE_COUNT + 3;

  /**
   * The feature id for the '<em><b>Target Class</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__TARGET_CLASS = MD_BUNDLE_MEMBER_FEATURE_COUNT + 4;

  /**
   * The feature id for the '<em><b>Description</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__DESCRIPTION = MD_BUNDLE_MEMBER_FEATURE_COUNT + 5;

  /**
   * The feature id for the '<em><b>Category</b></em>' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__CATEGORY = MD_BUNDLE_MEMBER_FEATURE_COUNT + 6;

  /**
   * The feature id for the '<em><b>Preview Image</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__PREVIEW_IMAGE = MD_BUNDLE_MEMBER_FEATURE_COUNT + 7;

  /**
   * The feature id for the '<em><b>Supported Features</b></em>' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__SUPPORTED_FEATURES = MD_BUNDLE_MEMBER_FEATURE_COUNT + 8;

  /**
   * The feature id for the '<em><b>Validator</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__VALIDATOR = MD_BUNDLE_MEMBER_FEATURE_COUNT + 9;

  /**
   * The feature id for the '<em><b>Supported Options</b></em>' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM__SUPPORTED_OPTIONS = MD_BUNDLE_MEMBER_FEATURE_COUNT + 10;

  /**
   * The number of structural features of the '<em>Md Algorithm</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_ALGORITHM_FEATURE_COUNT = MD_BUNDLE_MEMBER_FEATURE_COUNT + 11;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdCategoryImpl <em>Md Category</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdCategoryImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdCategory()
   * @generated
   */
  int MD_CATEGORY = 8;

  /**
   * The feature id for the '<em><b>Name</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY__NAME = MD_BUNDLE_MEMBER__NAME;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY__DOCUMENTATION = MD_BUNDLE_MEMBER__DOCUMENTATION;

  /**
   * The feature id for the '<em><b>Deprecated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY__DEPRECATED = MD_BUNDLE_MEMBER_FEATURE_COUNT + 0;

  /**
   * The feature id for the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY__LABEL = MD_BUNDLE_MEMBER_FEATURE_COUNT + 1;

  /**
   * The feature id for the '<em><b>Description</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY__DESCRIPTION = MD_BUNDLE_MEMBER_FEATURE_COUNT + 2;

  /**
   * The number of structural features of the '<em>Md Category</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_CATEGORY_FEATURE_COUNT = MD_BUNDLE_MEMBER_FEATURE_COUNT + 3;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionSupportImpl <em>Md Option Support</em>}' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionSupportImpl
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionSupport()
   * @generated
   */
  int MD_OPTION_SUPPORT = 9;

  /**
   * The feature id for the '<em><b>Option</b></em>' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_SUPPORT__OPTION = 0;

  /**
   * The feature id for the '<em><b>Value</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_SUPPORT__VALUE = 1;

  /**
   * The feature id for the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_SUPPORT__DOCUMENTATION = 2;

  /**
   * The number of structural features of the '<em>Md Option Support</em>' class.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   * @ordered
   */
  int MD_OPTION_SUPPORT_FEATURE_COUNT = 3;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.MdOptionTargetType <em>Md Option Target Type</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.MdOptionTargetType
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionTargetType()
   * @generated
   */
  int MD_OPTION_TARGET_TYPE = 10;

  /**
   * The meta object id for the '{@link org.eclipse.elk.core.meta.metaData.MdGraphFeature <em>Md Graph Feature</em>}' enum.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see org.eclipse.elk.core.meta.metaData.MdGraphFeature
   * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGraphFeature()
   * @generated
   */
  int MD_GRAPH_FEATURE = 11;


  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdModel <em>Md Model</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Model</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdModel
   * @generated
   */
  EClass getMdModel();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdModel#getName <em>Name</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Name</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdModel#getName()
   * @see #getMdModel()
   * @generated
   */
  EAttribute getMdModel_Name();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdModel#getImportSection <em>Import Section</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Import Section</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdModel#getImportSection()
   * @see #getMdModel()
   * @generated
   */
  EReference getMdModel_ImportSection();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdModel#getBundle <em>Bundle</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Bundle</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdModel#getBundle()
   * @see #getMdModel()
   * @generated
   */
  EReference getMdModel_Bundle();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdBundle <em>Md Bundle</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Bundle</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle
   * @generated
   */
  EClass getMdBundle();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getLabel <em>Label</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Label</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle#getLabel()
   * @see #getMdBundle()
   * @generated
   */
  EAttribute getMdBundle_Label();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getTargetClass <em>Target Class</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Target Class</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle#getTargetClass()
   * @see #getMdBundle()
   * @generated
   */
  EAttribute getMdBundle_TargetClass();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getDocumentationFolder <em>Documentation Folder</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Documentation Folder</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle#getDocumentationFolder()
   * @see #getMdBundle()
   * @generated
   */
  EAttribute getMdBundle_DocumentationFolder();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getIdPrefix <em>Id Prefix</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Id Prefix</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle#getIdPrefix()
   * @see #getMdBundle()
   * @generated
   */
  EAttribute getMdBundle_IdPrefix();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getMembers <em>Members</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Members</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle#getMembers()
   * @see #getMdBundle()
   * @generated
   */
  EReference getMdBundle_Members();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdBundleMember <em>Md Bundle Member</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Bundle Member</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundleMember
   * @generated
   */
  EClass getMdBundleMember();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundleMember#getName <em>Name</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Name</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundleMember#getName()
   * @see #getMdBundleMember()
   * @generated
   */
  EAttribute getMdBundleMember_Name();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdBundleMember#getDocumentation <em>Documentation</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Documentation</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdBundleMember#getDocumentation()
   * @see #getMdBundleMember()
   * @generated
   */
  EAttribute getMdBundleMember_Documentation();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdGroupOrOption <em>Md Group Or Option</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Group Or Option</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdGroupOrOption
   * @generated
   */
  EClass getMdGroupOrOption();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdGroup <em>Md Group</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Group</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdGroup
   * @generated
   */
  EClass getMdGroup();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.meta.metaData.MdGroup#getChildren <em>Children</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Children</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdGroup#getChildren()
   * @see #getMdGroup()
   * @generated
   */
  EReference getMdGroup_Children();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdOption <em>Md Option</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Option</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption
   * @generated
   */
  EClass getMdOption();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#isDeprecated <em>Deprecated</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Deprecated</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#isDeprecated()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Deprecated();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#isAdvanced <em>Advanced</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Advanced</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#isAdvanced()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Advanced();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#isProgrammatic <em>Programmatic</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Programmatic</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#isProgrammatic()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Programmatic();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#isOutput <em>Output</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Output</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#isOutput()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Output();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#isGlobal <em>Global</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Global</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#isGlobal()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Global();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOption#getType <em>Type</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Type</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getType()
   * @see #getMdOption()
   * @generated
   */
  EReference getMdOption_Type();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#getLabel <em>Label</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Label</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getLabel()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Label();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOption#getDescription <em>Description</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Description</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getDescription()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Description();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOption#getDefaultValue <em>Default Value</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Default Value</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getDefaultValue()
   * @see #getMdOption()
   * @generated
   */
  EReference getMdOption_DefaultValue();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOption#getLowerBound <em>Lower Bound</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Lower Bound</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getLowerBound()
   * @see #getMdOption()
   * @generated
   */
  EReference getMdOption_LowerBound();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOption#getUpperBound <em>Upper Bound</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Upper Bound</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getUpperBound()
   * @see #getMdOption()
   * @generated
   */
  EReference getMdOption_UpperBound();

  /**
   * Returns the meta object for the attribute list '{@link org.eclipse.elk.core.meta.metaData.MdOption#getTargets <em>Targets</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute list '<em>Targets</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getTargets()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_Targets();

  /**
   * Returns the meta object for the attribute list '{@link org.eclipse.elk.core.meta.metaData.MdOption#getLegacyIds <em>Legacy Ids</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute list '<em>Legacy Ids</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getLegacyIds()
   * @see #getMdOption()
   * @generated
   */
  EAttribute getMdOption_LegacyIds();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.meta.metaData.MdOption#getDependencies <em>Dependencies</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Dependencies</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOption#getDependencies()
   * @see #getMdOption()
   * @generated
   */
  EReference getMdOption_Dependencies();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdOptionDependency <em>Md Option Dependency</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Option Dependency</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionDependency
   * @generated
   */
  EClass getMdOptionDependency();

  /**
   * Returns the meta object for the reference '{@link org.eclipse.elk.core.meta.metaData.MdOptionDependency#getTarget <em>Target</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the reference '<em>Target</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionDependency#getTarget()
   * @see #getMdOptionDependency()
   * @generated
   */
  EReference getMdOptionDependency_Target();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOptionDependency#getValue <em>Value</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Value</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionDependency#getValue()
   * @see #getMdOptionDependency()
   * @generated
   */
  EReference getMdOptionDependency_Value();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm <em>Md Algorithm</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Algorithm</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm
   * @generated
   */
  EClass getMdAlgorithm();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#isDeprecated <em>Deprecated</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Deprecated</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#isDeprecated()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_Deprecated();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getProvider <em>Provider</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Provider</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getProvider()
   * @see #getMdAlgorithm()
   * @generated
   */
  EReference getMdAlgorithm_Provider();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getParameter <em>Parameter</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Parameter</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getParameter()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_Parameter();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getLabel <em>Label</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Label</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getLabel()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_Label();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getTargetClass <em>Target Class</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Target Class</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getTargetClass()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_TargetClass();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getDescription <em>Description</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Description</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getDescription()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_Description();

  /**
   * Returns the meta object for the reference '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getCategory <em>Category</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the reference '<em>Category</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getCategory()
   * @see #getMdAlgorithm()
   * @generated
   */
  EReference getMdAlgorithm_Category();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getPreviewImage <em>Preview Image</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Preview Image</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getPreviewImage()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_PreviewImage();

  /**
   * Returns the meta object for the attribute list '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedFeatures <em>Supported Features</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute list '<em>Supported Features</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedFeatures()
   * @see #getMdAlgorithm()
   * @generated
   */
  EAttribute getMdAlgorithm_SupportedFeatures();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getValidator <em>Validator</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Validator</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getValidator()
   * @see #getMdAlgorithm()
   * @generated
   */
  EReference getMdAlgorithm_Validator();

  /**
   * Returns the meta object for the containment reference list '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedOptions <em>Supported Options</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference list '<em>Supported Options</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedOptions()
   * @see #getMdAlgorithm()
   * @generated
   */
  EReference getMdAlgorithm_SupportedOptions();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdCategory <em>Md Category</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Category</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdCategory
   * @generated
   */
  EClass getMdCategory();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdCategory#isDeprecated <em>Deprecated</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Deprecated</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdCategory#isDeprecated()
   * @see #getMdCategory()
   * @generated
   */
  EAttribute getMdCategory_Deprecated();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdCategory#getLabel <em>Label</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Label</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdCategory#getLabel()
   * @see #getMdCategory()
   * @generated
   */
  EAttribute getMdCategory_Label();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdCategory#getDescription <em>Description</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Description</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdCategory#getDescription()
   * @see #getMdCategory()
   * @generated
   */
  EAttribute getMdCategory_Description();

  /**
   * Returns the meta object for class '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport <em>Md Option Support</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for class '<em>Md Option Support</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionSupport
   * @generated
   */
  EClass getMdOptionSupport();

  /**
   * Returns the meta object for the reference '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getOption <em>Option</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the reference '<em>Option</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionSupport#getOption()
   * @see #getMdOptionSupport()
   * @generated
   */
  EReference getMdOptionSupport_Option();

  /**
   * Returns the meta object for the containment reference '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getValue <em>Value</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the containment reference '<em>Value</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionSupport#getValue()
   * @see #getMdOptionSupport()
   * @generated
   */
  EReference getMdOptionSupport_Value();

  /**
   * Returns the meta object for the attribute '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getDocumentation <em>Documentation</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for the attribute '<em>Documentation</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionSupport#getDocumentation()
   * @see #getMdOptionSupport()
   * @generated
   */
  EAttribute getMdOptionSupport_Documentation();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.meta.metaData.MdOptionTargetType <em>Md Option Target Type</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Md Option Target Type</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionTargetType
   * @generated
   */
  EEnum getMdOptionTargetType();

  /**
   * Returns the meta object for enum '{@link org.eclipse.elk.core.meta.metaData.MdGraphFeature <em>Md Graph Feature</em>}'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the meta object for enum '<em>Md Graph Feature</em>'.
   * @see org.eclipse.elk.core.meta.metaData.MdGraphFeature
   * @generated
   */
  EEnum getMdGraphFeature();

  /**
   * Returns the factory that creates the instances of the model.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the factory that creates the instances of the model.
   * @generated
   */
  MetaDataFactory getMetaDataFactory();

  /**
   * <!-- begin-user-doc -->
   * Defines literals for the meta objects that represent
   * <ul>
   *   <li>each class,</li>
   *   <li>each feature of each class,</li>
   *   <li>each enum,</li>
   *   <li>and each data type</li>
   * </ul>
   * <!-- end-user-doc -->
   * @generated
   */
  interface Literals
  {
    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdModelImpl <em>Md Model</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdModelImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdModel()
     * @generated
     */
    EClass MD_MODEL = eINSTANCE.getMdModel();

    /**
     * The meta object literal for the '<em><b>Name</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_MODEL__NAME = eINSTANCE.getMdModel_Name();

    /**
     * The meta object literal for the '<em><b>Import Section</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_MODEL__IMPORT_SECTION = eINSTANCE.getMdModel_ImportSection();

    /**
     * The meta object literal for the '<em><b>Bundle</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_MODEL__BUNDLE = eINSTANCE.getMdModel_Bundle();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl <em>Md Bundle</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdBundleImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdBundle()
     * @generated
     */
    EClass MD_BUNDLE = eINSTANCE.getMdBundle();

    /**
     * The meta object literal for the '<em><b>Label</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE__LABEL = eINSTANCE.getMdBundle_Label();

    /**
     * The meta object literal for the '<em><b>Target Class</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE__TARGET_CLASS = eINSTANCE.getMdBundle_TargetClass();

    /**
     * The meta object literal for the '<em><b>Documentation Folder</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE__DOCUMENTATION_FOLDER = eINSTANCE.getMdBundle_DocumentationFolder();

    /**
     * The meta object literal for the '<em><b>Id Prefix</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE__ID_PREFIX = eINSTANCE.getMdBundle_IdPrefix();

    /**
     * The meta object literal for the '<em><b>Members</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_BUNDLE__MEMBERS = eINSTANCE.getMdBundle_Members();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdBundleMemberImpl <em>Md Bundle Member</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdBundleMemberImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdBundleMember()
     * @generated
     */
    EClass MD_BUNDLE_MEMBER = eINSTANCE.getMdBundleMember();

    /**
     * The meta object literal for the '<em><b>Name</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE_MEMBER__NAME = eINSTANCE.getMdBundleMember_Name();

    /**
     * The meta object literal for the '<em><b>Documentation</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_BUNDLE_MEMBER__DOCUMENTATION = eINSTANCE.getMdBundleMember_Documentation();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdGroupOrOptionImpl <em>Md Group Or Option</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdGroupOrOptionImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGroupOrOption()
     * @generated
     */
    EClass MD_GROUP_OR_OPTION = eINSTANCE.getMdGroupOrOption();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdGroupImpl <em>Md Group</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdGroupImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGroup()
     * @generated
     */
    EClass MD_GROUP = eINSTANCE.getMdGroup();

    /**
     * The meta object literal for the '<em><b>Children</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_GROUP__CHILDREN = eINSTANCE.getMdGroup_Children();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl <em>Md Option</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOption()
     * @generated
     */
    EClass MD_OPTION = eINSTANCE.getMdOption();

    /**
     * The meta object literal for the '<em><b>Deprecated</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__DEPRECATED = eINSTANCE.getMdOption_Deprecated();

    /**
     * The meta object literal for the '<em><b>Advanced</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__ADVANCED = eINSTANCE.getMdOption_Advanced();

    /**
     * The meta object literal for the '<em><b>Programmatic</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__PROGRAMMATIC = eINSTANCE.getMdOption_Programmatic();

    /**
     * The meta object literal for the '<em><b>Output</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__OUTPUT = eINSTANCE.getMdOption_Output();

    /**
     * The meta object literal for the '<em><b>Global</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__GLOBAL = eINSTANCE.getMdOption_Global();

    /**
     * The meta object literal for the '<em><b>Type</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION__TYPE = eINSTANCE.getMdOption_Type();

    /**
     * The meta object literal for the '<em><b>Label</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__LABEL = eINSTANCE.getMdOption_Label();

    /**
     * The meta object literal for the '<em><b>Description</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__DESCRIPTION = eINSTANCE.getMdOption_Description();

    /**
     * The meta object literal for the '<em><b>Default Value</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION__DEFAULT_VALUE = eINSTANCE.getMdOption_DefaultValue();

    /**
     * The meta object literal for the '<em><b>Lower Bound</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION__LOWER_BOUND = eINSTANCE.getMdOption_LowerBound();

    /**
     * The meta object literal for the '<em><b>Upper Bound</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION__UPPER_BOUND = eINSTANCE.getMdOption_UpperBound();

    /**
     * The meta object literal for the '<em><b>Targets</b></em>' attribute list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__TARGETS = eINSTANCE.getMdOption_Targets();

    /**
     * The meta object literal for the '<em><b>Legacy Ids</b></em>' attribute list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION__LEGACY_IDS = eINSTANCE.getMdOption_LegacyIds();

    /**
     * The meta object literal for the '<em><b>Dependencies</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION__DEPENDENCIES = eINSTANCE.getMdOption_Dependencies();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionDependencyImpl <em>Md Option Dependency</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionDependencyImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionDependency()
     * @generated
     */
    EClass MD_OPTION_DEPENDENCY = eINSTANCE.getMdOptionDependency();

    /**
     * The meta object literal for the '<em><b>Target</b></em>' reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION_DEPENDENCY__TARGET = eINSTANCE.getMdOptionDependency_Target();

    /**
     * The meta object literal for the '<em><b>Value</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION_DEPENDENCY__VALUE = eINSTANCE.getMdOptionDependency_Value();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl <em>Md Algorithm</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdAlgorithmImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdAlgorithm()
     * @generated
     */
    EClass MD_ALGORITHM = eINSTANCE.getMdAlgorithm();

    /**
     * The meta object literal for the '<em><b>Deprecated</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__DEPRECATED = eINSTANCE.getMdAlgorithm_Deprecated();

    /**
     * The meta object literal for the '<em><b>Provider</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_ALGORITHM__PROVIDER = eINSTANCE.getMdAlgorithm_Provider();

    /**
     * The meta object literal for the '<em><b>Parameter</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__PARAMETER = eINSTANCE.getMdAlgorithm_Parameter();

    /**
     * The meta object literal for the '<em><b>Label</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__LABEL = eINSTANCE.getMdAlgorithm_Label();

    /**
     * The meta object literal for the '<em><b>Target Class</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__TARGET_CLASS = eINSTANCE.getMdAlgorithm_TargetClass();

    /**
     * The meta object literal for the '<em><b>Description</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__DESCRIPTION = eINSTANCE.getMdAlgorithm_Description();

    /**
     * The meta object literal for the '<em><b>Category</b></em>' reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_ALGORITHM__CATEGORY = eINSTANCE.getMdAlgorithm_Category();

    /**
     * The meta object literal for the '<em><b>Preview Image</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__PREVIEW_IMAGE = eINSTANCE.getMdAlgorithm_PreviewImage();

    /**
     * The meta object literal for the '<em><b>Supported Features</b></em>' attribute list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_ALGORITHM__SUPPORTED_FEATURES = eINSTANCE.getMdAlgorithm_SupportedFeatures();

    /**
     * The meta object literal for the '<em><b>Validator</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_ALGORITHM__VALIDATOR = eINSTANCE.getMdAlgorithm_Validator();

    /**
     * The meta object literal for the '<em><b>Supported Options</b></em>' containment reference list feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_ALGORITHM__SUPPORTED_OPTIONS = eINSTANCE.getMdAlgorithm_SupportedOptions();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdCategoryImpl <em>Md Category</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdCategoryImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdCategory()
     * @generated
     */
    EClass MD_CATEGORY = eINSTANCE.getMdCategory();

    /**
     * The meta object literal for the '<em><b>Deprecated</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_CATEGORY__DEPRECATED = eINSTANCE.getMdCategory_Deprecated();

    /**
     * The meta object literal for the '<em><b>Label</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_CATEGORY__LABEL = eINSTANCE.getMdCategory_Label();

    /**
     * The meta object literal for the '<em><b>Description</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_CATEGORY__DESCRIPTION = eINSTANCE.getMdCategory_Description();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionSupportImpl <em>Md Option Support</em>}' class.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.impl.MdOptionSupportImpl
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionSupport()
     * @generated
     */
    EClass MD_OPTION_SUPPORT = eINSTANCE.getMdOptionSupport();

    /**
     * The meta object literal for the '<em><b>Option</b></em>' reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION_SUPPORT__OPTION = eINSTANCE.getMdOptionSupport_Option();

    /**
     * The meta object literal for the '<em><b>Value</b></em>' containment reference feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EReference MD_OPTION_SUPPORT__VALUE = eINSTANCE.getMdOptionSupport_Value();

    /**
     * The meta object literal for the '<em><b>Documentation</b></em>' attribute feature.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    EAttribute MD_OPTION_SUPPORT__DOCUMENTATION = eINSTANCE.getMdOptionSupport_Documentation();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.MdOptionTargetType <em>Md Option Target Type</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.MdOptionTargetType
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdOptionTargetType()
     * @generated
     */
    EEnum MD_OPTION_TARGET_TYPE = eINSTANCE.getMdOptionTargetType();

    /**
     * The meta object literal for the '{@link org.eclipse.elk.core.meta.metaData.MdGraphFeature <em>Md Graph Feature</em>}' enum.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.elk.core.meta.metaData.MdGraphFeature
     * @see org.eclipse.elk.core.meta.metaData.impl.MetaDataPackageImpl#getMdGraphFeature()
     * @generated
     */
    EEnum MD_GRAPH_FEATURE = eINSTANCE.getMdGraphFeature();

  }

} //MetaDataPackage
