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
package org.eclipse.elk.core.meta.metaData.impl;

import org.eclipse.elk.core.meta.metaData.*;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EDataType;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;

import org.eclipse.emf.ecore.impl.EFactoryImpl;

import org.eclipse.emf.ecore.plugin.EcorePlugin;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model <b>Factory</b>.
 * <!-- end-user-doc -->
 * @generated
 */
public class MetaDataFactoryImpl extends EFactoryImpl implements MetaDataFactory
{
  /**
   * Creates the default factory implementation.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public static MetaDataFactory init()
  {
    try
    {
      MetaDataFactory theMetaDataFactory = (MetaDataFactory)EPackage.Registry.INSTANCE.getEFactory(MetaDataPackage.eNS_URI);
      if (theMetaDataFactory != null)
      {
        return theMetaDataFactory;
      }
    }
    catch (Exception exception)
    {
      EcorePlugin.INSTANCE.log(exception);
    }
    return new MetaDataFactoryImpl();
  }

  /**
   * Creates an instance of the factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MetaDataFactoryImpl()
  {
    super();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public EObject create(EClass eClass)
  {
    switch (eClass.getClassifierID())
    {
      case MetaDataPackage.MD_MODEL: return createMdModel();
      case MetaDataPackage.MD_BUNDLE: return createMdBundle();
      case MetaDataPackage.MD_BUNDLE_MEMBER: return createMdBundleMember();
      case MetaDataPackage.MD_GROUP_OR_OPTION: return createMdGroupOrOption();
      case MetaDataPackage.MD_GROUP: return createMdGroup();
      case MetaDataPackage.MD_OPTION: return createMdOption();
      case MetaDataPackage.MD_OPTION_DEPENDENCY: return createMdOptionDependency();
      case MetaDataPackage.MD_ALGORITHM: return createMdAlgorithm();
      case MetaDataPackage.MD_CATEGORY: return createMdCategory();
      case MetaDataPackage.MD_OPTION_SUPPORT: return createMdOptionSupport();
      default:
        throw new IllegalArgumentException("The class '" + eClass.getName() + "' is not a valid classifier");
    }
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public Object createFromString(EDataType eDataType, String initialValue)
  {
    switch (eDataType.getClassifierID())
    {
      case MetaDataPackage.MD_OPTION_TARGET_TYPE:
        return createMdOptionTargetTypeFromString(eDataType, initialValue);
      case MetaDataPackage.MD_GRAPH_FEATURE:
        return createMdGraphFeatureFromString(eDataType, initialValue);
      default:
        throw new IllegalArgumentException("The datatype '" + eDataType.getName() + "' is not a valid classifier");
    }
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String convertToString(EDataType eDataType, Object instanceValue)
  {
    switch (eDataType.getClassifierID())
    {
      case MetaDataPackage.MD_OPTION_TARGET_TYPE:
        return convertMdOptionTargetTypeToString(eDataType, instanceValue);
      case MetaDataPackage.MD_GRAPH_FEATURE:
        return convertMdGraphFeatureToString(eDataType, instanceValue);
      default:
        throw new IllegalArgumentException("The datatype '" + eDataType.getName() + "' is not a valid classifier");
    }
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdModel createMdModel()
  {
    MdModelImpl mdModel = new MdModelImpl();
    return mdModel;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdBundle createMdBundle()
  {
    MdBundleImpl mdBundle = new MdBundleImpl();
    return mdBundle;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdBundleMember createMdBundleMember()
  {
    MdBundleMemberImpl mdBundleMember = new MdBundleMemberImpl();
    return mdBundleMember;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdGroupOrOption createMdGroupOrOption()
  {
    MdGroupOrOptionImpl mdGroupOrOption = new MdGroupOrOptionImpl();
    return mdGroupOrOption;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdGroup createMdGroup()
  {
    MdGroupImpl mdGroup = new MdGroupImpl();
    return mdGroup;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdOption createMdOption()
  {
    MdOptionImpl mdOption = new MdOptionImpl();
    return mdOption;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdOptionDependency createMdOptionDependency()
  {
    MdOptionDependencyImpl mdOptionDependency = new MdOptionDependencyImpl();
    return mdOptionDependency;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdAlgorithm createMdAlgorithm()
  {
    MdAlgorithmImpl mdAlgorithm = new MdAlgorithmImpl();
    return mdAlgorithm;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdCategory createMdCategory()
  {
    MdCategoryImpl mdCategory = new MdCategoryImpl();
    return mdCategory;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdOptionSupport createMdOptionSupport()
  {
    MdOptionSupportImpl mdOptionSupport = new MdOptionSupportImpl();
    return mdOptionSupport;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdOptionTargetType createMdOptionTargetTypeFromString(EDataType eDataType, String initialValue)
  {
    MdOptionTargetType result = MdOptionTargetType.get(initialValue);
    if (result == null) throw new IllegalArgumentException("The value '" + initialValue + "' is not a valid enumerator of '" + eDataType.getName() + "'");
    return result;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String convertMdOptionTargetTypeToString(EDataType eDataType, Object instanceValue)
  {
    return instanceValue == null ? null : instanceValue.toString();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MdGraphFeature createMdGraphFeatureFromString(EDataType eDataType, String initialValue)
  {
    MdGraphFeature result = MdGraphFeature.get(initialValue);
    if (result == null) throw new IllegalArgumentException("The value '" + initialValue + "' is not a valid enumerator of '" + eDataType.getName() + "'");
    return result;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String convertMdGraphFeatureToString(EDataType eDataType, Object instanceValue)
  {
    return instanceValue == null ? null : instanceValue.toString();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MetaDataPackage getMetaDataPackage()
  {
    return (MetaDataPackage)getEPackage();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @deprecated
   * @generated
   */
  @Deprecated
  public static MetaDataPackage getPackage()
  {
    return MetaDataPackage.eINSTANCE;
  }

} //MetaDataFactoryImpl
