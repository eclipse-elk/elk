/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
package org.eclipse.elk.core.debug.grandom.gRandom.util;

import org.eclipse.elk.core.debug.grandom.gRandom.*;

import org.eclipse.emf.common.notify.Adapter;
import org.eclipse.emf.common.notify.Notifier;

import org.eclipse.emf.common.notify.impl.AdapterFactoryImpl;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * The <b>Adapter Factory</b> for the model.
 * It provides an adapter <code>createXXX</code> method for each class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage
 * @generated
 */
public class GRandomAdapterFactory extends AdapterFactoryImpl
{
  /**
   * The cached model package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected static GRandomPackage modelPackage;

  /**
   * Creates an instance of the adapter factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public GRandomAdapterFactory()
  {
    if (modelPackage == null)
    {
      modelPackage = GRandomPackage.eINSTANCE;
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
  protected GRandomSwitch<Adapter> modelSwitch =
    new GRandomSwitch<Adapter>()
    {
      @Override
      public Adapter caseRandGraph(RandGraph object)
      {
        return createRandGraphAdapter();
      }
      @Override
      public Adapter caseConfiguration(Configuration object)
      {
        return createConfigurationAdapter();
      }
      @Override
      public Adapter caseHierarchy(Hierarchy object)
      {
        return createHierarchyAdapter();
      }
      @Override
      public Adapter caseEdges(Edges object)
      {
        return createEdgesAdapter();
      }
      @Override
      public Adapter caseNodes(Nodes object)
      {
        return createNodesAdapter();
      }
      @Override
      public Adapter caseSize(Size object)
      {
        return createSizeAdapter();
      }
      @Override
      public Adapter casePorts(Ports object)
      {
        return createPortsAdapter();
      }
      @Override
      public Adapter caseFlow(Flow object)
      {
        return createFlowAdapter();
      }
      @Override
      public Adapter caseDoubleQuantity(DoubleQuantity object)
      {
        return createDoubleQuantityAdapter();
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
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.RandGraph <em>Rand Graph</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.RandGraph
   * @generated
   */
  public Adapter createRandGraphAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration <em>Configuration</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Configuration
   * @generated
   */
  public Adapter createConfigurationAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy <em>Hierarchy</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy
   * @generated
   */
  public Adapter createHierarchyAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges <em>Edges</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Edges
   * @generated
   */
  public Adapter createEdgesAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Nodes <em>Nodes</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Nodes
   * @generated
   */
  public Adapter createNodesAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size <em>Size</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Size
   * @generated
   */
  public Adapter createSizeAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports <em>Ports</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Ports
   * @generated
   */
  public Adapter createPortsAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow <em>Flow</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Flow
   * @generated
   */
  public Adapter createFlowAdapter()
  {
    return null;
  }

  /**
   * Creates a new adapter for an object of class '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity <em>Double Quantity</em>}'.
   * <!-- begin-user-doc -->
   * This default implementation returns null so that we can easily ignore cases;
   * it's useful to ignore a case when inheritance will catch all the cases anyway.
   * <!-- end-user-doc -->
   * @return the new adapter.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity
   * @generated
   */
  public Adapter createDoubleQuantityAdapter()
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

} //GRandomAdapterFactory
