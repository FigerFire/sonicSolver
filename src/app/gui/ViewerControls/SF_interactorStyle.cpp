/// @file SF_interactorStyle.cpp
/// @brief 三维视图交互与工具栏控制实现。

#include "SF_interactorStyle.h"

#ifdef SONIC_GUI_WITH_VTK

#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

namespace SF::GUI {
namespace {

/// @brief ParaView 常用鼠标习惯：左键旋转，右键平移。
class ParaViewInteractorStyle final
    : public vtkInteractorStyleTrackballCamera {
public:
    static ParaViewInteractorStyle* New();
    vtkTypeMacro(ParaViewInteractorStyle, vtkInteractorStyleTrackballCamera);

    void OnLeftButtonDown() override
    {
        if (!Interactor) return;
        FindPokedRenderer(Interactor->GetEventPosition()[0],
                          Interactor->GetEventPosition()[1]);
        if (!CurrentRenderer) return;
        leftButtonDown_ = true;
        GrabFocus(EventCallbackCommand);
        StartRotate();
    }

    void OnLeftButtonUp() override
    {
        leftButtonDown_ = false;
        if (State == VTKIS_ROTATE) EndRotate();
        if (Interactor) ReleaseFocus();
        clearIdleState();
    }

    void OnMiddleButtonDown() override
    {
        if (!Interactor) return;
        FindPokedRenderer(Interactor->GetEventPosition()[0],
                          Interactor->GetEventPosition()[1]);
        if (!CurrentRenderer) return;
        middleButtonDown_ = true;
        GrabFocus(EventCallbackCommand);
        StartPan();
    }

    void OnMiddleButtonUp() override
    {
        middleButtonDown_ = false;
        if (State == VTKIS_PAN) EndPan();
        if (Interactor) ReleaseFocus();
        clearIdleState();
    }

    void OnRightButtonDown() override
    {
        if (!Interactor) return;
        FindPokedRenderer(Interactor->GetEventPosition()[0],
                          Interactor->GetEventPosition()[1]);
        if (!CurrentRenderer) return;
        rightButtonDown_ = true;
        GrabFocus(EventCallbackCommand);
        StartPan();
    }

    void OnRightButtonUp() override
    {
        rightButtonDown_ = false;
        if (State == VTKIS_PAN) EndPan();
        if (Interactor) ReleaseFocus();
        clearIdleState();
    }

    void OnMouseMove() override
    {
        // 普通鼠标划过 3D 窗口时不应该改变相机。
        // 这里只显式执行当前按键对应的动作，不交给基类猜内部状态。
        if (!anyButtonDown()) {
            clearIdleState();
            return;
        }
        if (!CurrentRenderer) {
            clearIdleState();
            return;
        }
        if (leftButtonDown_) {
            if (State != VTKIS_ROTATE) StartRotate();
            Rotate();
        } else if (rightButtonDown_ || middleButtonDown_) {
            if (State != VTKIS_PAN) StartPan();
            Pan();
        }
        // 每次交互都刷新裁剪范围和完整双缓冲帧，避免薄网格旋转后沿用旧深度。
        CurrentRenderer->ResetCameraClippingRange();
        if (Interactor) Interactor->Render();
    }

    void OnEnter() override
    {
        clearButtons();
        clearIdleState();
    }

    void OnLeave() override
    {
        clearButtons();
        clearIdleState();
    }

private:
    bool anyButtonDown() const
    {
        return leftButtonDown_ || middleButtonDown_ || rightButtonDown_;
    }

    void clearIdleState()
    {
        if (!anyButtonDown()) State = VTKIS_NONE;
    }

    void clearButtons()
    {
        leftButtonDown_ = false;
        middleButtonDown_ = false;
        rightButtonDown_ = false;
    }

    bool leftButtonDown_ = false;
    bool middleButtonDown_ = false;
    bool rightButtonDown_ = false;
};

vtkStandardNewMacro(ParaViewInteractorStyle);

} // namespace

vtkSmartPointer<vtkInteractorStyle> createParaViewInteractorStyle()
{
    return vtkSmartPointer<ParaViewInteractorStyle>::New();
}

} // namespace SF::GUI

#endif
