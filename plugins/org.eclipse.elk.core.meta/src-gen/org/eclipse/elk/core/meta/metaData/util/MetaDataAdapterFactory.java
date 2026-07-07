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
package org.eclipse.elk.core.meta.metaData.util;

import org.eclipse.elk.core.meta.metaData.*;

import org.eclipse.emf.common.notify.Adapter;
import org.eclipse.emf.common.notify.Notifier;

import org.eclipse.emf.common.notify.impl.AdapterFactoryImpl;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * The <b>Adapter Factory</b> for the model.
 * It provides an adapter <code>createXXX</code> method for each class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage
 * @generated
 */
public class MetaDataAdapterFactory extends AdapterFactoryImpl
{
  /**
   * The cached model package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected static MetaDataPackage modelPackage;

  /**
   * Creates an instance of the adapter factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MetaDataAdapterFactory()
  {
    if (modelPackage == null)
    {
      modelPackage = MetaDataPackage.eINSTANCE;
    }
  }

  /**
   * Returns whether this factory is applicable for the type of the object.
   * <!-- begin-user-doc -->
   * This implementation returns <code>true</code> if the object is either the model's package or is an instance object of the model.
   * <!-- end-user-doc -->
   * @return whether this factory is applicable for the type of the object.
   * @generated
   */
  @Override
  public boolean isFactoryForType(Object object)
  {
    if (object == modelPackage)
    {
      return true;
    }
    if (object instanceof EObject)
    {
      return ((EObject)object).eClass().getEPackage() == modelPackage;
    }
    return false;
  }

  /**
   * The switch that delegates to the <code>createXXX</code> methods.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected MetaDataSwitch<Adapter> modelSwitch =
    new MetaDataSwitch<Adapter>()
    {
      @Override
      public Adapter caseMdModel(MdModel object)
      {
        return createMdModelAdapter();
      }
      @Override
      public Adapter caseMdBundle(MdBundle object)
      {
        return createMdBundleAdapter();
      }
      @Override
      public Adapter caseMdBundleMember(MdBundleMember object)
      {
        return createMdBundleMemberAdapter();
      }
      @Override
      public Adapter caseMdGroupOrOption(MdGroupOrOption object)
      {
        return createMdGroupOrOptionAdapter();
      }
      @Override
      public Adapter caseMdGroup(MdGroup object)
      {
        return createMdGroupAdapter();
      }
      @Override
      public Adapter caseMdOption(MdOption object)
      {
        return createMdOptionAdapter();
      }
      @Override
      public Adapter caseMdOptionDependency(MdOptionDependency object)
      {
        return createMdOptionDependencyAdapter();
      }
      @Override
      public Adapter caseMdAlgorithm(MdAlgorithm object)
      {
        return createMdAlgorithmAdapter();
      }
      @Override
      public Adapter caseMdCategory(MdCategory object)
      {
        return createMdCategoryAdapter();
      }
      @Override
      public Adapter caseMdOptionSupport(MdOptionSupport object)
      {
        return createMdOptionSupportAdapter();
      }
      @Override
      public Adapter defaultCase(EObject object)
      {
        return createEObjectAdapter();
      }
    };

  /**
   * Creates an adapter for the <code>target</code>.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param target the object to adapt.
   * @return the adapter for the <code>target</code>.
   * @generated
   */
  @Override
  public Adapter createAdapter(Notifier target)
  {
    return modelSwitch.doSwitch((EObject)target);
  }


  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdModel <em>Md Model</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdModel
   * @generated
   */
  public Adapter createMdModelAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdBundle <em>Md Bundle</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdBundle
   * @generated
   */
  public Adapter createMdBundleAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdBundleMember <em>Md Bundle Member</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdBundleMember
   * @generated
   */
  public Adapter createMdBundleMemberAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdGroupOrOption <em>Md Group Or Option</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdGroupOrOption
   * @generated
   */
  public Adapter createMdGroupOrOptionAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdGroup <em>Md Group</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdGroup
   * @generated
   */
  public Adapter createMdGroupAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdOption <em>Md Option</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdOption
   * @generated
   */
  public Adapter createMdOptionAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdOptionDependency <em>Md Option Dependency</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionDependency
   * @generated
   */
  public Adapter createMdOptionDependencyAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm <em>Md Algorithm</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdAlgorithm
   * @generated
   */
  public Adapter createMdAlgorithmAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdCategory <em>Md Category</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdCategory
   * @generated
   */
  public Adapter createMdCategoryAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport <em>Md Option Support</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.meta.metaData.MdOptionSupport
   * @generated
   */
  public Adapter createMdOptionSupportAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for the default case.
   * <!-- begin-user-doc -->
   * This default implementation returns null.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @generated
   */
  public Adapter createEObjectAdapter()
  {
    return null;
  }

} //MetaDataAdapterFactory
